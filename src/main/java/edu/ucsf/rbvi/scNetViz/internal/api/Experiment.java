package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.List;
import javax.swing.table.TableModel;

// TODO: move this to API?
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;

public interface Experiment {
	public Source getSource();
	public Matrix getMatrix();
	public List<Category> getCategories();
	public Metadata getMetadata();
	public void addCategory(Category category);
	public String getSpecies();

	public Category getDefaultCategory();

	public DifferentialExpression getDiffExp();
	public void setDiffExp(DifferentialExpression expr);

	// For efficiency purposes, sometimes implementations
	// of Experiment might want to provide their own
	// TableModel
	public TableModel getTableModel();
}
