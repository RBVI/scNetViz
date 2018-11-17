package edu.ucsf.rbvi.scNetViz.internal.model;

/**
 * Approach to calculating differential expression
 *
 * 1) Filter all genes to remove all genes that don't have at least a 15% dDR (difference in detection rate)
 * 
 * 2) Calculate the differential expression between each gene in each cluster and every other cluster.  
 *    This is done by testing for differential gene expression for all genes above the dDR threshold in 
 *    every combination of clusters, then finding genes that have a positive gene expression ratio and 
 *    pass FDR threshold (default FDR < 1%) for a cluster in every comparison.
 * 3) Calculate the differential expression between genes in each cluster and /all/ other clusters.  This is
 *    the value that will be reported as the logGER.
 **/

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;

public class DifferentialExpression extends SimpleMatrix implements DoubleMatrix {
	final Category category;
	final Experiment experiment;
	final Matrix matrix;

	int nGenes;
	int nCategories;

	public DifferentialExpression(final ScNVManager manager, final Category category, int categoryRow) {
		super(manager);
		this.category = category;
		this.experiment = category.getExperiment();
		this.matrix = experiment.getMatrix();
		nGenes = matrix.getNRows();

	}

	// Calculate a list of genes (row numbers) to be ignored for all of the subsequent calculations
	// We do it as a list to avoid having to create a duplicate matrix
	public void filter(double dDRThreshold) {

	}

	public void findMarkers(double fdrThreshold) {
	}

	public void calculateDiffExp() {
	}

	@Override
	public String getMatrixType() { return "Simple String Matrix";}

	@Override
	public double getDoubleValue(int row, int column) {
		// return clusters[row][column];
		return 0.0;
	}

	@Override
	public double getDoubleValue(String row, String column) {
		// int col = colLabels.indexOf(column);

		// int intRow = Double.valueOf(row);
		// return clusters[intRow-minK][col-2];
		return 0.0;
	}

	@Override
	public double[][] getDoubleMatrix(double missingValue) {
		return null;
	}

}
