package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.Map;
import javax.swing.table.TableModel;

/*
 * A Category is a matrix of values organized with
 * the category key in the first column, and the
 * categorical values in each subsequent column.
 */
public interface Category {
	public String getCategoryType();
	public Matrix getMatrix();
	public Experiment getExperiment();
	public TableModel getTableModel();
	public int getHeaderCols();

	public double[][] getMeans(); // Where [][] = [nGenes][nCategories]
	public int[] getSizes(); 			// Where [] = [nCategories] and the contents are the number of cells in each category

	// dDRthreshold is the cutoff for the minimum difference between clusters
	public void filter(double dDRthreshold);

	// Calculate the logGER between each category and all other categories
	// This will trigger the calculation of means and sizes
	public Map<String, double[]> getLogGER();

	// Calculate the logGER between the category and all other categories
	// This will trigger the calculation of means and sizes
	public double[] getLogGER(String category1);

	// Calculate the logGER between the two categories
	// This will trigger the calculation of means and sizes
	public double[] getLogGER(String category1, String category2);

	public String toString();
}
