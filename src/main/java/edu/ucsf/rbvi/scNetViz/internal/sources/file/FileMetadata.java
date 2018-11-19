package edu.ucsf.rbvi.scNetViz.internal.sources.file;

import java.io.File;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;

public class FileMetadata extends HashMap<String, Object> implements Metadata {
	public static String ASSAYS = "assays";
	public static String FACTORS = "factors";
	public static String FILE = "file";

	public FileMetadata(File file) {
		super();
		put(ACCESSION, file.getName());
		put(TYPE,"File");
		put(DESCRIPTION, "Experiment loaded from file");
		put(FILE, file);
	}

	public String toString() {
		return "<html><p style='width: 500px'><b>"+get(ACCESSION)+"</b>: "+get(DESCRIPTION)+"</p></html>";
	}
}
