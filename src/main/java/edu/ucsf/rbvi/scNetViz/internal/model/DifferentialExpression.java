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
import java.util.List;
import java.util.Map;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class DifferentialExpression extends SimpleMatrix implements DoubleMatrix {
	final Category category;
	final Experiment experiment;
	int categoryRow;
	double dDRCutoff;
	double log2FCCutoff;
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
		Map<Object, Map<String, double[]>> logGERMap = category.getLogGER(categoryRow, dDRCutoff-.001, log2FCCutoff);
		super.nCols = means.keySet().size()*6; // 6 columns for each category/cluster

		setRowLabels(experiment.getMatrix().getRowLabels());

		// Get the column headers
		List<String> colHeaders = new ArrayList<>();
		colHeaders.add("Gene");
		for (Object cat: means.keySet()) {
			colHeaders.add(category.mkLabel(cat)+" MTC");
			colHeaders.add(category.mkLabel(cat)+" DR");
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
			// System.out.println("Means for "+cat+" col="+col);
			double[] mean = means.get(cat);
			double[] drs = drMap.get(cat);
			double[] mtdc = mtdcMap.get(cat);
			double[] logGER = logGERMap.get(cat).get("logFC");
			double[] pValue = logGERMap.get(cat).get("pValue");
			double[] FDR = adjustPValues(pValue);

			for (int row = 0; row < nRows; row++) {
				// System.out.println("row: "+row+" = "+mean[row]);
				matrix[col][row] = mean[row];
				matrix[col+1][row] = drs[row];
				matrix[col+2][row] = mtdc[row];
				matrix[col+3][row] = logGER[row];
				matrix[col+4][row] = pValue[row];
				matrix[col+5][row] = FDR[row];
			}
			// Arrays.fill(matrix[col+5], Double.NaN);
			col += 6;
		}
	}

	public SortableTableModel getTableModel() {
		if (tableModel == null)
			tableModel = new DiffExpTableModel(this, category, categoryRow);
		return tableModel;
	}

	// Calculate a list of genes (row numbers) to be ignored for all of the subsequent calculations
	// We do it as a list to avoid having to create a duplicate matrix
	public void filter(double dDRThreshold) {

	}

	public void findMarkers(double fdrThreshold) {
	}

	public void calculateDiffExp() {
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

	// Calculate the FDR using Benjamini Hochberg
	private double[] adjustPValues(double[] pValues) {
		double[] adjustedPvalues = new double[pValues.length];
		Arrays.fill(adjustedPvalues, Double.NaN);
		int testCount = countValues(pValues);
		Integer[] sortIndex = indexSort(pValues, pValues.length);
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

	public Integer[] indexSort(double[] tData, int nVals) {
    Integer[] index = new Integer[nVals];
    for (int i = 0; i < nVals; i++) index[i] = i;
    IndexComparator iCompare = new IndexComparator(tData);
    Arrays.sort(index, iCompare);
    return index;
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
		return getDoubleValue(intRow, col);
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
			return Double.class;
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

	  private static class IndexComparator implements Comparator<Integer> {
    double[] data = null;
    int[] intData = null;
    String[] stringData = null;

    public IndexComparator(String[] data) { this.stringData = data; }

    public IndexComparator(double[] data) { this.data = data; }

    public IndexComparator(int[] data) { this.intData = data; }

    public int compare(Integer o1, Integer o2) {
      if (data != null) {
				// NaN handling
				if (Double.isNaN(data[o1]) && Double.isNaN(data[o2]))
					return 0;
				if (Double.isNaN(data[o1]))
					return -1;
				if (Double.isNaN(data[o2]))
					return 1;

        if (data[o1] < data[o2]) return -1;
        if (data[o1] > data[o2]) return 1;
        return 0;
      } else if (intData != null) {
        if (intData[o1] < intData[o2]) return -1;
        if (intData[o1] > intData[o2]) return 1;
        return 0;
      } else if (stringData != null) {
        return stringData[o1].compareTo(stringData[o2]);
      }
      return 0;
    }
  }

}
