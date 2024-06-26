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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.json.simple.JSONObject;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.MyDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.PercentDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.PValueDouble;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVWriter;
import edu.ucsf.rbvi.scNetViz.internal.utils.MatrixUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class DifferentialExpression extends SimpleMatrix implements DoubleMatrix {
	final Category category;
	final Experiment experiment;
	int categoryRow;
	double dDRCutoff;
	double log2FCCutoff;
	Map<Object, Map<String, double[]>> logGERMap = null;
	// Map<Object, double[]> fdrMap = null;
	final double[][] matrix;

  List<String> colLabelCache = null;
  List<String> rowLabelCache = null;

	SortableTableModel tableModel = null;

	int nGenes;
	int nCategories;

	public DifferentialExpression(final ScNVManager manager, final Experiment experiment, 
	                              JSONObject json, File file) {
		super(manager);
		this.experiment = experiment;

		String catName = (String)json.get("category");
		Category catTemp = null;
		for (Category cat: experiment.getCategories()) {
			if (cat.toString().equals(catName)) {
				catTemp = cat;
				break;
			}
		}
		this.category = catTemp;

		List<String[]> matrixLines = null;
		try {
			matrixLines = CSVReader.readCSV(null, file);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		super.nCols = matrixLines.get(0).length-1;
		super.setColLabels(Arrays.asList(matrixLines.get(0)),0);
		super.nRows = matrixLines.size()-1;
		matrix = new double[nRows][nCols];
		List<String> rowLabels = new ArrayList<>();
		for (int row = 1; row < matrixLines.size(); row++) {
			String[] line = matrixLines.get(row);
			rowLabels.add(line[0]);
			for (int col = 1; col < line.length; col++) {
				try {
					if (line[col].length() == 0)
						matrix[row-1][col-1] = Double.NaN;
					else
						matrix[row-1][col-1] = Double.parseDouble(line[col]);
				} catch (NumberFormatException nfe) {
					throw new RuntimeException("Unable to read line: "+line[col]);
				}
			}
		}
		super.setRowLabels(rowLabels, 0);

		// Now, built our logGERMap
		logGERMap = new HashMap<>();
		for (Object cat: category.getMeans(categoryRow).keySet()) {
			double[] logGer = new double[nRows];
			double[] pValue = new double[nRows];
			double[] FDR = new double[nRows];

			int gerColumn = getColumn(cat, "log2FC");
			int pVColumn = getColumn(cat, "pValue");
			int fdrColumn = getColumn(cat, "FDR");

			// System.out.println("gerColumn = "+gerColumn+", pVColumn = "+pVColumn+", fdrColumn = "+fdrColumn);

			for (int row = 0; row < nRows; row++) {
				// System.out.println("log2FC = "+getDoubleValue(row, gerColumn)+", pValue = "+getDoubleValue(row, pVColumn));
				logGer[row] = getDoubleValue(row, gerColumn-1);
				pValue[row] = getDoubleValue(row, pVColumn-1);
				FDR[row] = getDoubleValue(row, fdrColumn-1);
			}
			
			HashMap<String, double[]> logCatMap = new HashMap<>();
			logCatMap.put("logFC", logGer);
			logCatMap.put("pValue", pValue);
			logCatMap.put("FDR", FDR);
			logGERMap.put(cat, logCatMap);

		}
	}

	public DifferentialExpression(final ScNVManager manager, final Category category, int categoryRow,
	                              double dDRCutoff, double log2FCCutoff) {
		super(manager);
		this.category = category;
		this.categoryRow = categoryRow;
		this.dDRCutoff = dDRCutoff;
		this.log2FCCutoff = log2FCCutoff;
		this.experiment = category.getExperiment();
		Matrix mtx = experiment.getMatrix();
		super.nRows = experiment.getMatrix().getNRows();
		if (mtx instanceof MatrixMarket) {
			super.nRows = super.nRows - ((MatrixMarket)mtx).findControls().cardinality();
		}

		long start = System.currentTimeMillis();

		Map<Object, double[]> means = category.getMeans(categoryRow);
		System.out.println("getMeans took "+(System.currentTimeMillis()-start));
		long mid = System.currentTimeMillis();
		Map<Object, double[]> drMap = category.getDr(categoryRow);
		System.out.println("getDr took "+(System.currentTimeMillis()-mid));
		mid = System.currentTimeMillis();
		Map<Object, double[]> mtdcMap = category.getMTDC(categoryRow);
		System.out.println("mtdcMap took "+(System.currentTimeMillis()-mid));
		// fdrMap = new HashMap<>();

		mid = System.currentTimeMillis();
		logGERMap = category.getLogGER(categoryRow, dDRCutoff-.001, log2FCCutoff);
		System.out.println("getLogGER took "+(System.currentTimeMillis()-mid));
		super.nCols = means.keySet().size()*6; // 6 columns for each category/cluster
		if (means.containsKey(Category.UNUSED_CAT))
			super.nCols = super.nCols - 6;  // We don't want to show the "unused" category

		List<String> labels = new ArrayList<String>(nRows);
		for (int row = 0; row < mtx.getNRows(); row++) {
			if (mtx instanceof MatrixMarket && ((MatrixMarket)mtx).isControl(row))
				continue;
			labels.add(mtx.getRowLabels().get(row)[0]);
		}
		setRowLabels(labels, 0);

		// Get the column headers
		List<String> colHeaders = new ArrayList<>();
		colHeaders.add("Gene");
		for (Object cat: means.keySet()) {
			if (cat.equals(Category.UNUSED_CAT))
				continue;
			colHeaders.add(category.mkLabel(cat)+" MTC");
			colHeaders.add(category.mkLabel(cat)+" Min.pct");
			colHeaders.add(category.mkLabel(cat)+" MDTC");
			colHeaders.add(category.mkLabel(cat)+" log2FC");
			colHeaders.add(category.mkLabel(cat)+" pValue");
			colHeaders.add(category.mkLabel(cat)+" FDR");
		}
		setColLabels(colHeaders, 0);

		// Initialize the matrix
		matrix = new double[nRows][nCols];
		int col = 0;
		for (Object cat: means.keySet()) {
			if (cat.equals(Category.UNUSED_CAT))
				continue;

			double[] mean = means.get(cat);
			double[] drs = drMap.get(cat);
			double[] mtdc = mtdcMap.get(cat);
			double[] logGER = logGERMap.get(cat).get("logFC");
			double[] pValue = logGERMap.get(cat).get("pValue");
			double[] FDR = adjustPValues(pValue);
			logGERMap.get(cat).put("FDR", FDR);
			// fdrMap.put(cat, FDR);

			for (int row = 0; row < nRows; row++) {
				matrix[row][col] = mean[row];
				matrix[row][col+1] = drs[row];
				matrix[row][col+2] = mtdc[row];
				matrix[row][col+3] = logGER[row];
				matrix[row][col+4] = pValue[row];
				matrix[row][col+5] = FDR[row];
			}
			col += 6;
		}
		System.out.println("Total time: "+(System.currentTimeMillis()-start));
	}

	public Experiment getExperiment() { return experiment; }

	public int getCategoryRow() { return categoryRow; }
	public Category getCategory() { return category; }

	public SortableTableModel getTableModel() {
		if (tableModel == null) {
			tableModel = new DiffExpTableModel(this, category, categoryRow);
		}
		return tableModel;
	}

	public String toJSON() {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		builder.append("\"category\": \""+category.toString()+"\",");
		builder.append("\"row\": "+categoryRow+",");
		builder.append("\"ddrCutoff\": "+dDRCutoff+",");
		builder.append("\"log2FCCutoff\": "+log2FCCutoff+"}");
		return builder.toString();
	}

  public List<String> getRowLabels(int lbl) {
    return super.getRowLabels(lbl);
  }

  public double getLog2FCCutoff() { return log2FCCutoff; }

	public String toString() {
		return "Differential expression for category "+category+" row "+categoryRow;
	}

	public Map<Object,Map<String, double[]>> getLogGERMap() { 
		return logGERMap; 
	}

	public Map<String, double[]> getLogGERMap(Object cat) { 
		return logGERMap.get(cat); 
	}

	public double[] getLogGER(Object cat, boolean positiveOnly) { 
		if (logGERMap.containsKey(cat)) {
			double[] allGER = logGERMap.get(cat).get("logFC"); 

			if (positiveOnly) {
				double[] posGER = Arrays.copyOf(allGER, allGER.length);
				for (int i = 0; i < posGER.length; i++) {
					if (!Double.isNaN(posGER[i]) && posGER[i] < 0.0)
						posGER[i] = Double.NaN;
				}
				return posGER;
			} else {
				return allGER;
			}
		}
		return null;
	}

	public double[] getLogGER(Object cat) { 
		if (logGERMap.containsKey(cat))
			return logGERMap.get(cat).get("logFC"); 
		return null;
	}

	public List<String> getMarkerGenes(Object cat, double ddRCutoff, double log2FCCutoff, int minCells) {
		return null;
	}

	public double[] getGeneList(Object cat, double fdrCutoff, double log2FCCutoff, int nGenes, 
                                  boolean positiveOnly, int maxGenes, final List<String> geneList) {
		double[] logGER = logGERMap.get(cat).get("logFC");
		double[] fdr = logGERMap.get(cat).get("FDR");
		//
		// System.out.println("getGeneList: nGenes = "+nGenes+", maxGenes = "+maxGenes+", positiveOnly = "+positiveOnly);

		// if (nGenes > 0)
		// 	return getTopGenes(null, logGER, nGenes);

		double fc[] = new double[nRows];
		Arrays.fill(fc, Double.NaN);
		int count = 0;

		for (int row = 0; row < nRows; row++) {
			double thisfc = logGER[row];
			if (!positiveOnly)
				thisfc = Math.abs(thisfc);

			if (thisfc > log2FCCutoff && fdr[row] < fdrCutoff) {
				geneList.add(getRowLabel(row));
				fc[count++] = logGER[row];
			}
		}
		if (count > maxGenes) {
			return getTopGenes(geneList, fc, maxGenes);
		} else {
			return getTopGenes(geneList, fc, count);
		}
	}

	public Category getCurrentCategory() {
		return category;
	}

	public Set<Object> getCategoryValues() {
		// Make sure we exclude our "unused" category
		if (category.getMeans(categoryRow).containsKey(Category.UNUSED_CAT)) {
			Set<Object> set = new HashSet<>(category.getMeans(categoryRow).keySet());
			set.remove(Category.UNUSED_CAT);
			return set;
		}
		return category.getMeans(categoryRow).keySet();
	}

	public void saveFile(File file) throws IOException {
		CSVWriter.writeCSV(file, this, "\t");
	}

	private int getColumn(Object cat, String lbl) {
		String colName = category.mkLabel(cat)+" "+lbl;
		for (int col = 1; col <= nCols; col++) {
			if (colName.equals(getColumnLabel(col)))
				return col;
		}
		return -1;
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

	private double[] getTopGenes(List<String> geneList, double[] fc, int nGenes) {
		// Sort
		Integer[] sortedFC = MatrixUtils.indexSort(fc, fc.length, true);

		/*
		// Skip over the NaN's
		int start = 0;
		for (start = 0; start < fc.length; start++) {
			if (!Double.isNaN(fc[sortedFC[start]]))
				break;
		}

		List<String> newGeneList = new ArrayList<String>();
		if (geneList == null) {
			for (int topGene = 0; topGene < nGenes; topGene++) {
				newGeneList.add(getRowLabel(sortedFC[topGene+start]));
			}
		} else {
			for (int topGene = 0; topGene < nGenes; topGene++) {
				newGeneList.add(geneList.get(sortedFC[topGene+start]));
			}
		}
		*/
		double[] returnFC = new double[nGenes];
		List<String> newGeneList = new ArrayList<String>();
		for (int topGene = 0; topGene < nGenes; topGene++) {
			newGeneList.add(geneList.get(sortedFC[sortedFC.length-topGene-1]));
			returnFC[topGene] = fc[sortedFC[sortedFC.length-topGene-1]];
		}
		geneList.clear();
		geneList.addAll(newGeneList);
		return returnFC;
	}

	// Calculate the FDR using Benjamini Hochberg
	private double[] adjustPValues(double[] pValues) {
		double[] adjustedPvalues = new double[pValues.length];
		Arrays.fill(adjustedPvalues, Double.NaN);
		int testCount = countValues(pValues);
		// System.out.println("  testCount = "+testCount);
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
		return matrix[row][column];
	}

	@Override
	public double getDoubleValue(String row, String column) {
    if (colLabelCache == null) {
      colLabelCache = getColLabels(columnKey);
    }
		int col = colLabelCache.indexOf(column);
    if (rowLabelCache == null) {
      rowLabelCache = getRowLabels(rowKey);
    }
		int intRow = rowLabelCache.indexOf(row);
		if (col < 0 || intRow < 0) {
			System.out.println("Got negative index for "+row+","+column);
		}
		return getDoubleValue(intRow, col-1);
	}

	@Override
	public double[][] getDoubleMatrix(double missingValue) {
		return getDoubleMatrix(missingValue, false);
	}

	@Override
	public double[][] getDoubleMatrix(double missingValue, boolean transpose) {
		return getDoubleMatrix(missingValue, transpose, true);
	}

	// Note: we've already excluded controls
	@Override
	public double[][] getDoubleMatrix(double missingValue, boolean transpose, boolean excludeControls) {
		if (transpose && transposed)
			transpose = false;
		if (!transpose)
			return matrix;

		double[][] newArray;
	 	if (transpose)
			newArray = new double[nCols][nRows];
		else
			newArray = new double[nRows][nCols];

		for (int row = 0; row < nRows; row++) {
			for (int col = 0; col < nCols; col++) {
				if (transpose)
					newArray[col][row] = matrix[row][col];
				else
					newArray[row][col] = matrix[row][col];
			}
		}
		return newArray;
	}


  @Override
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
