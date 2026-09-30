package org.biojava.tornado;

import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.CifFileReader;

/** Test structures bundled as resources, so the tests need no network. */
public final class TestStructures {

	private TestStructures() {
	}

	public static Structure load(String id) throws Exception {
		try (InputStream in = new GZIPInputStream(
				TestStructures.class.getResourceAsStream("/pdb/" + id.toLowerCase() + ".cif.gz"))) {
			return new CifFileReader().getStructure(in);
		}
	}

	public static Atom[] caAtoms(String id, String chain) throws Exception {
		return StructureTools.getRepresentativeAtomArray(load(id).getPolyChainByPDB(chain));
	}

	/** Runs the given code with the GPU path forced on (or, without a TornadoVM runtime, the fallback). */
	public static <T> T forced(ThrowingSupplier<T> code) throws Exception {
		System.setProperty(TornadoSupport.PROPERTY, "force");
		try {
			return code.get();
		} finally {
			System.clearProperty(TornadoSupport.PROPERTY);
		}
	}

	public interface ThrowingSupplier<T> {
		T get() throws Exception;
	}
}
