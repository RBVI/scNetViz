package edu.ucsf.rbvi.scNetViz.internal.api;

import java.io.File;
import java.util.List;
import java.util.Map;
import javax.swing.table.TableModel;

// TODO: move this to API?
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;

public interface Experiment {
	public Source getSource();
	public Matrix getMatrix();
	public List<Category> getCategories();
	public Category getCategory(String categoryName);
	public Metadata getMetadata();
	public void addCategory(Category category);
	public String getSpecies();

	public Category getDefaultCategory();

	public DifferentialExpression getDiffExp();
	public void setDiffExp(DifferentialExpression expr);

	public String toJSON();
	public String toHTML();

	public Map<String,double[]> getTSNE();
	public void setTSNE(Map<String,double[]> tnse);
	public String getPlotType();
	public void setPlotType(String type);

	// For efficiency purposes, sometimes implementations
	// of Experiment might want to provide their own
	// TableModel
	public TableModel getTableModel();

	// This is the hook for saving all of the experiment files in a session
	public void createSessionFiles(String accession, List<File> files) throws Exception;
}
