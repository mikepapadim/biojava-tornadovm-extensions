package org.biojava.tornado.asa;

import uk.ac.manchester.tornado.api.KernelContext;
import uk.ac.manchester.tornado.api.types.arrays.DoubleArray;
import uk.ac.manchester.tornado.api.types.arrays.IntArray;

/**
 * GPU kernels for the Shrake-Rupley accessible surface area, following {@code AsaCalculator}.
 * <p>
 * Everything is computed in double precision, as in BioJava, so the per-atom results match exactly.
 */
public final class AsaKernels {

	/** Work-group size: one work-group per atom. */
	public static final int GROUP_SIZE = 128;

	private AsaKernels() {
	}

	/** Maximum number of neighbours per atom held in local memory by {@link #fusedAsa}. */
	public static final int LOCAL_CAPACITY = 160;

	/**
	 * Fused neighbour search, sort and occlusion test: work-group {@code g} computes the number of accessible
	 * sphere points of atom {@code g}, keeping its neighbour list in local memory.
	 * <ol>
	 * <li>The 27 grid cells around the atom are laid out as one virtual candidate list; in rounds of GROUP_SIZE,
	 * every thread tests one candidate with the {@code AsaCalculator.areNeighbors} criterion and the hits are
	 * compacted into local memory with a work-group prefix sum (deterministic, no atomics).</li>
	 * <li>The neighbours are rank-sorted by increasing distance (ties by position), as
	 * {@code AsaCalculator.Neighbors.createSorted}; the order only matters for the speed of step 3.</li>
	 * <li>Each thread tests sphere points against the sorted neighbours, stopping at the first occluding one, and
	 * the counts are summed in local memory.</li>
	 * </ol>
	 * If an atom has more than {@link #LOCAL_CAPACITY} neighbours its count is set to -1 and the host falls back to
	 * the CPU.
	 *
	 * @param ctx kernel context
	 * @param coords atom coordinates, interleaved x,y,z
	 * @param radii van der Waals radii
	 * @param cellStart start of each grid cell in cellAtoms, size nCells+1
	 * @param cellAtoms atom indices sorted by cell
	 * @param grid [minX, minY, minZ, cellSize, probe]
	 * @param dims [nAtoms, nx, ny, nz]
	 * @param spherePoints unit sphere points, interleaved x,y,z
	 * @param counts output, number of accessible points per atom, or -1 on neighbour overflow
	 * @param nPoints number of sphere points
	 */
	public static void fusedAsa(KernelContext ctx, DoubleArray coords, DoubleArray radii, IntArray cellStart,
			IntArray cellAtoms, DoubleArray grid, IntArray dims, DoubleArray spherePoints, IntArray counts,
			int nPoints) {
		int i = ctx.groupIdx;
		int lid = ctx.localIdx;
		int[] scan = ctx.allocateIntLocalArray(GROUP_SIZE);
		int[] cellBegin = ctx.allocateIntLocalArray(32);
		int[] cellPrefix = ctx.allocateIntLocalArray(32);
		double[] found = ctx.allocateDoubleLocalArray(4 * LOCAL_CAPACITY);
		double[] sorted = ctx.allocateDoubleLocalArray(4 * LOCAL_CAPACITY);
		double[] keys = ctx.allocateDoubleLocalArray(LOCAL_CAPACITY);

		int nx = dims.get(1);
		int ny = dims.get(2);
		int nz = dims.get(3);
		double cellSize = grid.get(3);
		double probe = grid.get(4);
		double xi = coords.get(3 * i);
		double yi = coords.get(3 * i + 1);
		double zi = coords.get(3 * i + 2);
		double vdwI = radii.get(i);
		// same operation order as AsaCalculator, so that the results match to the last bit
		double ri = probe + vdwI;

		// 1a. the 27 candidate cells as one virtual list
		if (lid < 27) {
			int gx = (int) ((xi - grid.get(0)) / cellSize) + lid / 9 - 1;
			int gy = (int) ((yi - grid.get(1)) / cellSize) + (lid / 3) % 3 - 1;
			int gz = (int) ((zi - grid.get(2)) / cellSize) + lid % 3 - 1;
			int begin = 0;
			int length = 0;
			if (gx >= 0 && gx < nx && gy >= 0 && gy < ny && gz >= 0 && gz < nz) {
				int cell = (gx * ny + gy) * nz + gz;
				begin = cellStart.get(cell);
				length = cellStart.get(cell + 1) - begin;
			}
			cellBegin[lid] = begin;
			cellPrefix[lid + 1] = length;
		}
		ctx.localBarrier();
		if (lid == 0) {
			cellPrefix[0] = 0;
			for (int c = 0; c < 27; c++) {
				cellPrefix[c + 1] += cellPrefix[c];
			}
		}
		ctx.localBarrier();
		int total = cellPrefix[27];

		// 1b. test candidates in rounds, compacting the hits with a prefix sum
		int count = 0;
		for (int roundBase = 0; roundBase < total; roundBase += GROUP_SIZE) {
			int c = roundBase + lid;
			int hit = 0;
			double dx = 0;
			double dy = 0;
			double dz = 0;
			double cutoff = 0;
			if (c < total) {
				int k = 0;
				while (cellPrefix[k + 1] <= c) {
					k++;
				}
				int j = cellAtoms.get(cellBegin[k] + c - cellPrefix[k]);
				if (j != i) {
					dx = coords.get(3 * j) - xi;
					dy = coords.get(3 * j + 1) - yi;
					dz = coords.get(3 * j + 2) - zi;
					double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
					// AsaCalculator.areNeighbors
					if (dist < vdwI + probe + probe + radii.get(j)) {
						hit = 1;
						double rj = radii.get(j) + probe;
						// equation 3 in Eisenhaber 1994
						cutoff = (dist * dist + ri * ri - rj * rj) / (2 * ri);
					}
				}
			}
			scan[lid] = hit;
			for (int offset = 1; offset < GROUP_SIZE; offset *= 2) {
				ctx.localBarrier();
				int add = lid >= offset ? scan[lid - offset] : 0;
				ctx.localBarrier();
				scan[lid] += add;
			}
			ctx.localBarrier();
			int slot = count + scan[lid] - 1;
			if (hit == 1 && slot < LOCAL_CAPACITY) {
				found[4 * slot] = dx;
				found[4 * slot + 1] = dy;
				found[4 * slot + 2] = dz;
				found[4 * slot + 3] = cutoff;
			}
			count += scan[GROUP_SIZE - 1];
			ctx.localBarrier();
		}

		if (count > LOCAL_CAPACITY) {
			if (lid == 0) {
				counts.set(i, -1);
			}
		} else {
			// 2. rank sort by squared distance, ties by position
			for (int k = lid; k < count; k += GROUP_SIZE) {
				double x = found[4 * k];
				double y = found[4 * k + 1];
				double z = found[4 * k + 2];
				keys[k] = x * x + y * y + z * z;
			}
			ctx.localBarrier();
			for (int k = lid; k < count; k += GROUP_SIZE) {
				double key = keys[k];
				int rank = 0;
				for (int m = 0; m < count; m++) {
					double other = keys[m];
					if (other < key || (other == key && m < k)) {
						rank++;
					}
				}
				sorted[4 * rank] = found[4 * k];
				sorted[4 * rank + 1] = found[4 * k + 1];
				sorted[4 * rank + 2] = found[4 * k + 2];
				sorted[4 * rank + 3] = found[4 * k + 3];
			}
			ctx.localBarrier();

			// 3. occlusion test, as AsaCalculator.calcSingleAsa
			int accessible = 0;
			for (int p = lid; p < nPoints; p += GROUP_SIZE) {
				double px = spherePoints.get(3 * p);
				double py = spherePoints.get(3 * p + 1);
				double pz = spherePoints.get(3 * p + 2);
				int occluded = 0;
				for (int k = 0; k < count && occluded == 0; k++) {
					double dotProd = sorted[4 * k] * px + sorted[4 * k + 1] * py + sorted[4 * k + 2] * pz;
					if (dotProd > sorted[4 * k + 3]) {
						occluded = 1;
					}
				}
				accessible += 1 - occluded;
			}
			scan[lid] = accessible;
			for (int stride = GROUP_SIZE / 2; stride > 0; stride /= 2) {
				ctx.localBarrier();
				if (lid < stride) {
					scan[lid] += scan[lid + stride];
				}
			}
			if (lid == 0) {
				counts.set(i, scan[0]);
			}
		}
	}
}
