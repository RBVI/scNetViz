package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.List;

public interface Experiment {
	public Source getSource();
	public Matrix getMatrix();
	public List<Category> getCategories();
	public Metadata getMetadata();
}
