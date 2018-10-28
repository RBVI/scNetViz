package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.List;
import javax.swing.table.TableModel;

public interface Experiment {
	public Source getSource();
	public Matrix getMatrix();
	public List<Category> getCategories();
	public Metadata getMetadata();

	// For efficiency purposes, sometimes implementations
	// of Experiment might want to provide their own
	// TableModel
	public TableModel getTableModel();
}
