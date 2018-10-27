package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.List;

public interface Source {
	public String getName();
	public List<String> getAccessions();
	public List<Metadata> getMetadata();
	public Experiment getExperiment(String accession);
}
