package edu.ucsf.rbvi.scNetViz.internal.api;

import java.io.File;
import java.util.List;
import java.util.Map;
import org.json.simple.JSONObject;

public interface Source {
	public static String SOURCENAME = "name";
	public String getName();
	public List<String> getAccessions();
	public List<Metadata> getMetadata();
	public Experiment getExperiment(String accession);
	public Experiment loadExperimentFromSession(JSONObject jsonExperiment, Map<String,File> fileMap);
	public Category loadCategoryFromSession(JSONObject jsonCategory, Experiment experiment, Map<String,File> fileMap);
}
