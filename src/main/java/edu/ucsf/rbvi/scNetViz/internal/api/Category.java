package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.Map;
import javax.swing.table.TableModel;

/*
 * A Category is a matrix of values organized with
 * the category key in the first column, and the
 * categorical values in each subsequent column.
 */
public interface Category {
	
	/**
	 * Return the category "type", which is actually just
	 * a string description of the category.
	 *
	 * @return category type
	 */
	public String getCategoryType();

	/**
	 * Get the Matrix for this category.
	 *
	 * @return category matrix
	 */
	public Matrix getMatrix();

	/**
	 * Return the Experiment that this category refers to
	 *
	 * @return the Experiment for this category
	 */
	public Experiment getExperiment();

	/**
	 * Return a table model suitable for visualizing this category
	 * in a JTable
	 *
	 * @return the table model for this category
	 */
	public TableModel getTableModel();

	/**
	 * Return the number of header columns.
	 *
	 * @return the number of header columns
	 */
	public int getHeaderCols();

	/**
	 * Get the currently selected row (if any)
	 *
	 * @return the row number of the selected row
	 */
	public int getSelectedRow();

	/**
	 * Set the currently selected row
	 *
	 * @param selectedRow
	 */
	public void setSelectedRow(int selectedRow);

	/**
	 * Return the mean values for every gene and each unique value
	 * for a category.
	 *
	 * @param category this is the category row we're using to get the means
	 * @return a Map with the category value as the key and a double array of
	 * means for each gene across that category
	 */
	public Map<Object, double[]> getMeans(int category);

	/**
	 * Return the dr values for every gene and each unique value
	 * for a category.
	 *
	 * @param category this is the category row we're using to get the means
	 * @return a Map with the category value as the key and a double array of
	 * the proportion of this gene in the category
	 */
	public Map<Object, double[]> getDr(int category);

	public Map<Object, double[]> getMTDC(int category);

	/**
	 * Return the sizes (number of assays) for each unique value
	 * for a category.
	 *
	 * @param category this is the category row we're using to get the sizes
	 * @return a Map with the category value as the key and the size as the value
	 */
	public Map<Object, Integer> getSizes(int category);

	// dDRthreshold is the cutoff for the minimum difference between clusters
	public void filter(int category, double dDRthreshold);

	// Calculate the logGER between each category and all other categories
	// This will trigger the calculation of means and sizes
	public Map<String, double[]> getLogGER(int category);

	// Calculate the logGER between the category and all other categories
	// This will trigger the calculation of means and sizes
	public double[] getLogGER(int category, String category1);

	// Calculate the logGER between the two categories
	// This will trigger the calculation of means and sizes
	public double[] getLogGER(int category, String category1, String category2);

	public String toString();
}
