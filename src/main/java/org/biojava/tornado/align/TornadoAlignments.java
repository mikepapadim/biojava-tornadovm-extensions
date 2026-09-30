package org.biojava.tornado.align;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.GuideTree;
import org.biojava.nbio.alignment.template.AbstractScorer;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.tornado.TornadoSupport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.manchester.tornado.api.GridScheduler;
import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.TaskGraph;
import uk.ac.manchester.tornado.api.TornadoExecutionPlan;
import uk.ac.manchester.tornado.api.WorkerGrid;
import uk.ac.manchester.tornado.api.WorkerGrid1D;
import uk.ac.manchester.tornado.api.enums.DataTransferMode;
import uk.ac.manchester.tornado.api.types.arrays.IntArray;

/**
 * GPU versions of the all-pairs scoring entry points of {@link Alignments}.
 * <p>
 * Supports the {@link PairwiseSequenceScorerType#GLOBAL GLOBAL} (Needleman-Wunsch) and
 * {@link PairwiseSequenceScorerType#LOCAL LOCAL} (Smith-Waterman) score types with affine or constant gap penalties;
 * the scores are exactly BioJava's. Other score types (identities and similarities need a traceback), linear gap
 * penalties, small inputs, and runs without a TornadoVM runtime use {@link Alignments}.
 */
public final class TornadoAlignments {

	private static final Logger logger = LoggerFactory.getLogger(TornadoAlignments.class);

	/** Below this many DP cells in total the CPU implementation is used. */
	public static final long DEFAULT_MIN_GPU_CELLS = 50_000_000L;

	/** Device memory budget for the DP row buffers of one batch, in bytes. */
	static final long BUFFER_BUDGET = 1L << 30;

	private TornadoAlignments() {
	}

	/** @return true if the given scorer type is computed on the GPU */
	public static boolean isSupported(PairwiseSequenceScorerType type, GapPenalty gapPenalty) {
		return (type == PairwiseSequenceScorerType.GLOBAL || type == PairwiseSequenceScorerType.LOCAL)
				&& gapPenalty.getType() != GapPenalty.Type.LINEAR;
	}

	/**
	 * Drop-in for {@link Alignments#getAllPairsScores(List, PairwiseSequenceScorerType, GapPenalty,
	 * SubstitutionMatrix)}: the score of every pair (i, j), i &lt; j, in the same order.
	 */
	public static <S extends Sequence<C>, C extends Compound> double[] getAllPairsScores(List<S> sequences,
			PairwiseSequenceScorerType type, GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		if (!isSupported(type, gapPenalty) || !TornadoSupport.useGpu(totalCells(sequences), DEFAULT_MIN_GPU_CELLS)) {
			return Alignments.getAllPairsScores(sequences, type, gapPenalty, subMatrix);
		}
		try {
			int[] scores = computeScores(sequences, type == PairwiseSequenceScorerType.LOCAL, gapPenalty, subMatrix);
			double[] out = new double[scores.length];
			for (int k = 0; k < scores.length; k++) {
				out[k] = scores[k];
			}
			return out;
		} catch (RuntimeException | Error e) {
			TornadoSupport.disable(e);
			return Alignments.getAllPairsScores(sequences, type, gapPenalty, subMatrix);
		}
	}

	/**
	 * Drop-in for {@link Alignments#getAllPairsScorers} followed by {@link Alignments#runPairwiseScorers}: scorers
	 * whose scores are already computed, with the same {@link PairwiseSequenceScorer#getMaxScore() max} and
	 * {@link PairwiseSequenceScorer#getMinScore() min} (hence distances) as BioJava's aligners. They can be passed
	 * to {@link GuideTree#GuideTree(List, List)}.
	 */
	public static <S extends Sequence<C>, C extends Compound> List<PairwiseSequenceScorer<S, C>> getAllPairsScorers(
			List<S> sequences, PairwiseSequenceScorerType type, GapPenalty gapPenalty,
			SubstitutionMatrix<C> subMatrix) {
		if (!isSupported(type, gapPenalty) || !TornadoSupport.useGpu(totalCells(sequences), DEFAULT_MIN_GPU_CELLS)) {
			List<PairwiseSequenceScorer<S, C>> scorers = Alignments.getAllPairsScorers(sequences, type, gapPenalty,
					subMatrix);
			Alignments.runPairwiseScorers(scorers);
			return scorers;
		}
		boolean local = type == PairwiseSequenceScorerType.LOCAL;
		int[] scores;
		try {
			scores = computeScores(sequences, local, gapPenalty, subMatrix);
		} catch (RuntimeException | Error e) {
			TornadoSupport.disable(e);
			return getAllPairsScorers(sequences, type, gapPenalty, subMatrix);
		}
		int[] selfScores = new int[sequences.size()];
		for (int i = 0; i < sequences.size(); i++) {
			for (C c : sequences.get(i)) {
				selfScores[i] += subMatrix.getValue(c, c);
			}
		}
		List<PairwiseSequenceScorer<S, C>> scorers = new ArrayList<>(scores.length);
		int k = 0;
		for (int i = 0; i < sequences.size(); i++) {
			for (int j = i + 1; j < sequences.size(); j++) {
				S query = sequences.get(i);
				S target = sequences.get(j);
				// as AbstractPairwiseSequenceAligner.reset
				int max = Math.max(selfScores[i], selfScores[j]);
				int min = local ? 0 : (int) (2 * gapPenalty.getOpenPenalty()
						+ (query.getLength() + target.getLength()) * gapPenalty.getExtensionPenalty());
				scorers.add(new PrecomputedScorer<>(query, target, scores[k++], max, min));
			}
		}
		return scorers;
	}

	/** A {@link PairwiseSequenceScorer} holding an already computed score. */
	public static final class PrecomputedScorer<S extends Sequence<C>, C extends Compound> extends AbstractScorer
			implements PairwiseSequenceScorer<S, C> {
		private final S query;
		private final S target;
		private final int score;
		private final int max;
		private final int min;

		PrecomputedScorer(S query, S target, int score, int max, int min) {
			this.query = query;
			this.target = target;
			this.score = score;
			this.max = max;
			this.min = min;
		}

		@Override
		public S getQuery() {
			return query;
		}

		@Override
		public S getTarget() {
			return target;
		}

		@Override
		public double getMaxScore() {
			return max;
		}

		@Override
		public double getMinScore() {
			return min;
		}

		@Override
		public double getScore() {
			return score;
		}
	}

	private static <S extends Sequence<C>, C extends Compound> long totalCells(List<S> sequences) {
		long sum = 0, sumSq = 0;
		for (S s : sequences) {
			sum += s.getLength();
			sumSq += (long) s.getLength() * s.getLength();
		}
		return (sum * sum - sumSq) / 2;
	}

	/**
	 * Encodes the sequences, orders the pairs by DP size (so that the threads of a warp do similar work) and runs
	 * them in batches that fit the buffer budget.
	 */
	static <S extends Sequence<C>, C extends Compound> int[] computeScores(List<S> sequences, boolean local,
			GapPenalty gapPenalty, SubstitutionMatrix<C> subMatrix) {
		int nSeq = sequences.size();
		int nPairs = nSeq * (nSeq - 1) / 2;
		int[] result = new int[nPairs];
		if (nPairs == 0) {
			return result;
		}

		// alphabet: every distinct compound, scored with the matrix's own getValue so that unknown compounds get
		// exactly the value BioJava gives them
		Map<C, Integer> codes = new HashMap<>();
		List<C> alphabet = new ArrayList<>();
		int[] seqStart = new int[nSeq + 1];
		int maxLength = 0;
		for (int i = 0; i < nSeq; i++) {
			int length = sequences.get(i).getLength();
			seqStart[i + 1] = seqStart[i] + length;
			maxLength = Math.max(maxLength, length);
		}
		int[] residues = new int[Math.max(1, seqStart[nSeq])];
		for (int i = 0; i < nSeq; i++) {
			int pos = seqStart[i];
			for (C c : sequences.get(i)) {
				Integer code = codes.get(c);
				if (code == null) {
					code = alphabet.size();
					codes.put(c, code);
					alphabet.add(c);
				}
				residues[pos++] = code;
			}
		}
		int k = alphabet.size();
		int[] subs = new int[Math.max(1, k * k)];
		for (int a = 0; a < k; a++) {
			for (int b = 0; b < k; b++) {
				subs[a * k + b] = subMatrix.getValue(alphabet.get(a), alphabet.get(b));
			}
		}

		// pairs in BioJava's order, then sorted by decreasing DP size
		long[] keyed = new long[nPairs];
		int[] pairQuery = new int[nPairs];
		int[] pairTarget = new int[nPairs];
		int idx = 0;
		for (int i = 0; i < nSeq; i++) {
			for (int j = i + 1; j < nSeq; j++) {
				pairQuery[idx] = i;
				pairTarget[idx] = j;
				long cells = (long) sequences.get(i).getLength() * sequences.get(j).getLength();
				// high bits: DP size (descending), low bits: pair index
				keyed[idx] = ((Long.MAX_VALUE >> 24) - Math.min(cells, Long.MAX_VALUE >> 24)) << 24 | idx;
				idx++;
			}
		}
		if (nPairs >= (1 << 24)) {
			// the pair index no longer fits the sort key: skip the ordering
			for (int p = 0; p < nPairs; p++) {
				keyed[p] = p;
			}
		} else {
			Arrays.sort(keyed);
		}
		int[] order = new int[nPairs];
		for (int p = 0; p < nPairs; p++) {
			order[p] = (int) (keyed[p] & ((1 << 24) - 1));
		}

		int rowBytes = 3 * 4 * (maxLength + 1);
		int batch = (int) Math.min(nPairs, Math.max(SwKernels.GROUP_SIZE, BUFFER_BUDGET / rowBytes));
		batch = (batch + SwKernels.GROUP_SIZE - 1) / SwKernels.GROUP_SIZE * SwKernels.GROUP_SIZE;

		Batcher batcher = new Batcher(residues, seqStart, subs, k, gapPenalty.getOpenPenalty(),
				gapPenalty.getExtensionPenalty(), local, batch, maxLength);
		try {
			for (int from = 0; from < nPairs; from += batch) {
				int size = Math.min(batch, nPairs - from);
				int[] scores = batcher.run(order, pairQuery, pairTarget, from, size);
				for (int p = 0; p < size; p++) {
					result[order[from + p]] = scores[p];
				}
			}
		} finally {
			batcher.close();
		}
		return result;
	}

	/** One execution plan, re-executed for every batch of pairs. */
	private static final class Batcher implements AutoCloseable {
		private final int batch;
		private final IntArray pairQuery;
		private final IntArray pairTarget;
		private final IntArray params;
		private final IntArray scores;
		private final WorkerGrid worker;
		private final GridScheduler gridScheduler;
		private final TornadoExecutionPlan plan;

		Batcher(int[] residues, int[] seqStart, int[] subs, int alphabet, int gop, int gep, boolean local,
				int batch, int maxLength) {
			this.batch = batch;
			IntArray residueArray = IntArray.fromArray(residues);
			IntArray seqStartArray = IntArray.fromArray(seqStart);
			IntArray subsArray = IntArray.fromArray(subs);
			pairQuery = new IntArray(batch);
			pairTarget = new IntArray(batch);
			params = IntArray.fromElements(alphabet, gop, gep, local ? 1 : 0, 0, batch);
			scores = new IntArray(batch);
			long bufferSize = (long) (maxLength + 1) * batch;
			IntArray bufM = new IntArray((int) bufferSize);
			IntArray bufS0 = new IntArray((int) bufferSize);
			IntArray bufS1 = new IntArray((int) bufferSize);
			KernelContext ctx = new KernelContext();

			TaskGraph taskGraph = new TaskGraph("sw")
					.transferToDevice(DataTransferMode.FIRST_EXECUTION, residueArray, seqStartArray, subsArray)
					.transferToDevice(DataTransferMode.EVERY_EXECUTION, pairQuery, pairTarget, params)
					.task("scores", SwKernels::alignmentScores, ctx, residueArray, seqStartArray, pairQuery,
							pairTarget, subsArray, params, bufM, bufS0, bufS1, scores)
					.transferToHost(DataTransferMode.EVERY_EXECUTION, scores);
			worker = new WorkerGrid1D(batch);
			worker.setLocalWork(SwKernels.GROUP_SIZE, 1, 1);
			gridScheduler = new GridScheduler("sw.scores", worker);
			plan = new TornadoExecutionPlan(taskGraph.snapshot());
		}

		int[] run(int[] order, int[] allQuery, int[] allTarget, int from, int size) {
			for (int p = 0; p < size; p++) {
				int pair = order[from + p];
				pairQuery.set(p, allQuery[pair]);
				pairTarget.set(p, allTarget[pair]);
			}
			params.set(SwKernels.P_PAIRS, size);
			worker.setGlobalWork((size + SwKernels.GROUP_SIZE - 1) / SwKernels.GROUP_SIZE * SwKernels.GROUP_SIZE,
					1, 1);
			long t0 = System.nanoTime();
			plan.withGridScheduler(gridScheduler).execute();
			if (Boolean.getBoolean("biojava.tornado.trace")) {
				logger.info("sw: batch of {} pairs in {} ms", size, (System.nanoTime() - t0) / 1e6);
			}
			int[] out = new int[size];
			for (int p = 0; p < size; p++) {
				out[p] = scores.get(p);
			}
			return out;
		}

		@Override
		public void close() {
			try {
				plan.close();
			} catch (Exception e) {
				logger.warn("Could not release the alignment execution plan", e);
			}
		}
	}
}
