package org.biojava.tornado.bench;

import java.io.File;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.io.CifFileReader;

/** Loads benchmark structures from {@code data/pdb/<id>.cif.gz} when present, else through {@link StructureIO}. */
final class Structures {

	private Structures() {
	}

	static Structure load(String id) throws Exception {
		File f = new File(System.getProperty("bench.pdbdir", "data/pdb"), id.toLowerCase() + ".cif.gz");
		if (f.isFile()) {
			return new CifFileReader().getStructure(f);
		}
		return StructureIO.getStructure(id);
	}
}
