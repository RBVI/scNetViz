package edu.ucsf.rbvi.scNetViz.internal.model;

/**
 * Approach to calculating differential expression
 *
 * 1) Filter all genes to remove all genes that don't have at least a 10% dDR (difference in detection rate)
 * 
 * 2) Calculate the differential expression between each gene in each cluster and every other cluster.  
 *    This is done by testing for differential gene expression for all genes above the dDR threshold in 
 *    every combination of clusters, then finding genes that have a positive gene expression ratio and 
 *    pass FDR threshold (default FDR < 1%) for a cluster in every comparison.
 * 3) Calculate the differential expression between genes in each cluster and /all/ other clusters.  This is
 *    the value that will be reported as the logGER.
 **/

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.MyDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.PercentDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.PValueDouble;
import edu.ucsf.rbvi.scNetViz.internal.utils.MatrixUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class DifferentialExpression extends SimpleMatrix implements DoubleMatrix {
	final Category category;
	final Experiment experiment;
	int categoryRow;
	double dDRCutoff;
	double log2FCCutoff;
	Map<Object, Map<String, double[]>> logGERMap;
	Map<Object, double[]> fdrMap;
	final double[][] matrix;

	SortableTableModel tableModel = null;

	int nGenes;
	int nCategories;

	public DifferentialExpression(final ScNVManager manager, final Category category, int categoryRow,
	                              double dDRCutoff, double log2FCCutoff) {
		super(manager);
		this.category = category;
		this.categoryRow = categoryRow;
		this.dDRCutoff = dDRCutoff;
		this.log2FCCutoff = log2FCCutoff;
		this.experiment = category.getExperiment();
		super.nRows = experiment.getMatrix().getNRows();
		Map<Object, double[]> means = category.getMeans(categoryRow);
		Map<Object, double[]> drMap = category.getDr(categoryRow);
		Map<Object, double[]> mtdcMap = category.getMTDC(categoryRow);
		fdrMap = new HashMap<>();

		logGERMap = category.getLogGER(categoryRow, dDRCutoff-.001, log2FCCutoff);
		super.nCols = means.keySet().size()*6; // 6 columns for each category/cluster

		setRowLabels(experiment.getMatrix().getRowLabels());

		// Get the column headers
		List<String> colHeaders = new ArrayList<>();
		colHeaders.add("Gene");
		for (Object cat: means.keySet()) {
			colHeaders.add(category.mkLabel(cat)+" MTC");
			colHeaders.add(category.mkLabel(cat)+" Min.pct");
			colHeaders.add(category.mkLabel(cat)+" MDTC");
			colHeaders.add(category.mkLabel(cat)+" logGER");
			colHeaders.add(category.mkLabel(cat)+" pValue");
			colHeaders.add(category.mkLabel(cat)+" FDR");
		}
		setColLabels(colHeaders);

		// Initialize the matrix
		matrix = new double[nCols][nRows];
		int col = 0;
		for (Object cat: means.keySet()) {
			double[] mean = means.get(cat);
			double[] drs = drMap.get(cat);
			double[] mtdc = mtdcMap.get(cat);
			double[] logGER = logGERMap.get(cat).get("logFC");
			double[] pValue = logGERMap.get(cat).get("pValue");
			double[] FDR = adjustPValues(pValue);
			fdrMap.put(cat, FDR);

			for (int row = 0; row < nRows; row++) {
				// System.out.println("row: "+row+" = "+mean[row]);
				matrix[col][row] = mean[row];
				matrix[col+1][row] = drs[row];
				matrix[col+2][row] = mtdc[row];
				matrix[col+3][row] = logGER[row];
				matrix[col+4][row] = pValue[row];
				matrix[col+5][row] = FDR[row];
			}
			col += 6;
		}
	}

	public SortableTableModel getTableModel() {
		if (tableModel == null)
			tableModel = new DiffExpTableModel(this, category, categoryRow);
		return tableModel;
	}

	public List<String> getGeneList(Object cat, double fdrCutoff, double log2FCCutoff, int nGenes, int maxGenes) {
		double[] logGER = logGERMap.get(cat).get("logFC");
		double[] pValues = logGERMap.get(cat).get("pValue");
		double[] fdr = fdrMap.get(cat);

		List<String> geneList = new ArrayList<>();
		if (nGenes > 0) {
			return getTopGenes(geneList, pValues, nGenes);
		} else {
			// Should this be pValue or fdr?
			double pV[] = new double[nRows];
			Arrays.fill(pV, Double.NaN);
			int count = 0;
			for (int row = 0; row < nRows; row++) {
				if (Math.abs(logGER[row]) > log2FCCutoff && fdr[row] < fdrCutoff) {
					geneList.add(getRowLabel(row));
					pV[count++] = pValues[row];
				}
			}
			if (count > maxGenes) {
				geneList.clear();
				return getTopGenes(geneList, pV, maxGenes);
			}
		}
		return geneList;
	}

	public Category getCurrentCategory() {
		return category;
	}

	public Set<Object> getCategoryValues() {
		return category.getMeans(categoryRow).keySet();
	}

	private int countValues(double[] array) {
		int count = 0;
		for (int index = 0; index < array.length; index++) {
			if (!Double.isNaN(array[index]))
				count++;
		}
		return count;
	}

	// Simple Bonferroni adjustment
	private double adjustPValue(double pValue, int testCount) {
		return Math.min(1.0, pValue*(double)testCount);
	}

	private List<String> getTopGenes(List<String> geneList, double[] pValues, int nGenes) {
		// Sort
		Integer[] sortedPValues = MatrixUtils.indexSort(pValues, pValues.length);

		// Skip over the NaN's
		int start = 0;
		for (start = 0; start < pValues.length; start++) {
			if (!Double.isNaN(pValues[sortedPValues[start]]))
				break;
		}

		for (int topGene = 0; topGene < nGenes; topGene++) {
			geneList.add(getRowLabel(sortedPValues[topGene+start]));
		}
		return geneList;
	}

	// Calculate the FDR using Benjamini Hochberg
	private double[] adjustPValues(double[] pValues) {
		double[] adjustedPvalues = new double[pValues.length];
		Arrays.fill(adjustedPvalues, Double.NaN);
		int testCount = countValues(pValues);
		System.out.println("  testCount = "+testCount);
		Integer[] sortIndex = MatrixUtils.indexSort(pValues, pValues.length);
		/*
		for (int index = (pValues.length-testCount); index < pValues.length; index++) {
			System.out.println("pValues["+sortIndex[index]+"] = "+pValues[sortIndex[index]]);
		}
		*/

		for (int i = pValues.length-1; i >=(pValues.length-testCount); i--) {
			int index = sortIndex[i];
			if (i == pValues.length-1) {
				adjustedPvalues[index] = pValues[index];
			} else {
				double unadjustedPvalue = pValues[index];
				int divideByM = i+1-pValues.length+testCount;
				double left = adjustedPvalues[sortIndex[i+1]];
				double right = (testCount / (double) divideByM) * unadjustedPvalue;
				adjustedPvalues[index] = Math.min(left, right);
			}
		}

		return adjustedPvalues;
	}

	@Override
	public String getMatrixType() { return "Simple String Matrix";}

	@Override
	public double getDoubleValue(int row, int column) {
		return matrix[column][row];
	}

	@Override
	public double getDoubleValue(String row, String column) {
		int col = colLabels.indexOf(column);
		int intRow = rowLabels.indexOf(row);
		return getDoubleValue(intRow, col-1);
	}

	@Override
	public double[][] getDoubleMatrix(double missingValue) {
		return matrix;
	}

	public Double getValue(int row, int column) {
		double v = getDoubleValue(row, column);
		if (Double.isNaN(v)) return null;
		return v;
	}

	public class DiffExpTableModel extends SortableTableModel {
		final DifferentialExpression diffExp;
		final Category category;
		final Experiment experiment;
		int categoryRow;
		String sortedRow;

		DiffExpTableModel(final DifferentialExpression diffexp,
		                  final Category category, int categoryRow) {
			super(1);
			this.category = category;
			this.experiment = category.getExperiment();
			this.categoryRow = categoryRow;
			this.diffExp = diffexp;
		}

		@Override
		public int getColumnCount() { return diffExp.getNCols()+1; }

		@Override
		public String getColumnName(int column) {
			if (columnIndex == null) 
				return strip(diffExp.getColumnLabel(column));
			else
				return strip(diffExp.getColumnLabel(columnIndex[column]));
		}

		@Override
		public int getRowCount() { 
			return diffExp.getNRows();
		}

		@Override
		public Class<?> getColumnClass(int column) {
			if (column == 0)
				return String.class;

			int col = column%6;
			switch(col) {
				case 1:
				case 3:
				case 4:
					return MyDouble.class;
				case 2:
					return PercentDouble.class;
				case 5:
				case 0:
					return PValueDouble.class;
				default:
					return Double.class;
			}
		}

		@Override
		public Object getValueAt(int row, int column) {
			if (column == 0) {
				return strip(diffExp.getRowLabel(row));
			}

			if (columnIndex != null)
				column = columnIndex[column];

			Double v = diffExp.getValue(row, column-1);
			// System.out.println("getValueAt("+row+","+column+") = "+v);
			return v;
		}

		@Override
		public void sortColumns(int row) {
			sortedRow = strip(diffExp.getRowLabel(row));
			super.sortColumns(row);
		}

		@Override
		public int getSelectedRow() { return category.getSelectedRow(); }

		@Override
		public void setSelectedRow(int selectedRow) { category.setSelectedRow(selectedRow); }

		public String strip(String str) {
			return str.replaceAll("^\"|\"$", "");
		}
	}

}
