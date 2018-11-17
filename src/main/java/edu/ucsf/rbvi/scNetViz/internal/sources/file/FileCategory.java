package edu.ucsf.rbvi.scNetViz.internal.sources.file;

import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import org.apache.log4j.Logger;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.StringMatrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.SimpleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.utils.LogUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class FileCategory extends SimpleMatrix implements Category, StringMatrix {
	final Logger logger;
	final Experiment experiment;
	final String name;
	int hdrCols = 1;

	String[][] categories;

	String sortedRow = null;

	SortableTableModel tableModel = null;

	public FileCategory(final ScNVManager scManager, final Experiment experiment, final String name) {
		super(scManager);
		this.experiment = experiment;
		this.name = name;
		logger = Logger.getLogger(CyUserLog.NAME);
	}

	@Override
	public String toString() { return name;}

	@Override
	public String getCategoryType() { return name;}

	@Override
	public Experiment getExperiment() { return experiment;}

	@Override
	public String getMatrixType() { return "Simple String Matrix";}

	@Override
	public Matrix getMatrix() { return this;}

	@Override
	public String[][] getStringMatrix() { return categories;}

	@Override
	public String getValue(int row, int col) { return categories[row][col];}

	@Override
	public String getValue(String rowLabel, String colLabel) { 
		int row = rowLabels.indexOf(rowLabel);
		int col = colLabels.indexOf(colLabel);
		return categories[row][col];
	}

	@Override
	public int getHeaderCols() { return 1; }

	@Override
	public double[][] getMeans() {
		// Where [][] = [nGenes][nCategories]
		return null;
	}

	@Override
	public int[] getSizes() {
		// Where [] = [nCategories] and the contents are the number of cells in each category
		return null;
	}

	// dDRthreshold is the cutoff for the minimum difference between clusters
	@Override
	public void filter(double dDRthreshold) {
		return;
	}

	// Calculate the logGER between each category and all other categories
	// This will trigger the calculation of means and sizes
	@Override
	public Map<String, double[]> getLogGER() {
		return null;
	}

	// Calculate the logGER between the category and all other categories
	// This will trigger the calculation of means and sizes
	@Override
	public double[] getLogGER(String category1) {
		return null;
	};

	// Calculate the logGER between the two categories
	// This will trigger the calculation of means and sizes
	@Override
	public double[] getLogGER(String category1, String category2) {
		return null;
	}

	public static FileCategory fetchCategory(ScNVManager scManager, Experiment experiment,
	                                         File file, boolean transpose, int hdrCols,
	                                         TaskMonitor monitor) throws Exception {

		System.out.println("fetchCategory: file = "+file.toString()+", transpose = "+transpose);

		List<String[]> input = CSVReader.readCSV(monitor, file);
		if (input == null || input.size() < 2) {
			System.out.println("No input!");
			return null;
		}

		FileCategory fileCategory = new FileCategory(scManager, experiment, file.getName());
		fileCategory.hdrCols = hdrCols;

		List<String> labels;

		if (!transpose) {
			fileCategory.nRows = input.size();
			fileCategory.nCols = input.get(0).length-1;
			fileCategory.setColLabels(Arrays.asList(input.get(0)));
			labels = new ArrayList<String>(fileCategory.nCols);
		} else {
			fileCategory.nCols = input.size();
			fileCategory.nRows = input.get(0).length-1;
			fileCategory.setRowLabels(Arrays.asList(input.get(0)));
			labels = new ArrayList<String>(fileCategory.nRows);
		}

		fileCategory.categories = new String[fileCategory.nCols][fileCategory.nRows];

		labels.add("Category");

		boolean first = true;
		int lineNumber = 0;
		for (String[] line: input) {
			if (first) {
				first = false;
			} else {
				labels.add(line[0]);
				if (!transpose) {
					for (int col = 1; col < fileCategory.nCols; col++) {
						fileCategory.categories[col][lineNumber] = line[col];
					}
				} else {
					for (int row = 1; row < fileCategory.nRows; row++) {
						fileCategory.categories[lineNumber][row] = line[row];
					}
				}
			}
			lineNumber++;
		}

		if (!transpose) {
			fileCategory.setRowLabels(labels);
		} else {
			fileCategory.setColLabels(labels);
		}

		LogUtils.log(monitor, TaskMonitor.Level.INFO, "Read "+fileCategory.nRows+
			                    " rows with "+fileCategory.nCols+" columns");

		System.out.println("Read "+fileCategory.nRows+" rows with "+fileCategory.nCols+" columns");
		return fileCategory;
	}

	public Map<String,List<String>> getClusterList(String factor) {
		Map<String, List<String>> clusterMap = new HashMap<>();

		int factorColumn = colLabels.indexOf(factor);

		for (int row = 0; row < nRows; row++) {
			String id = rowLabels.get(row);
			String rowFactor = categories[row][factorColumn];
			if (!clusterMap.containsKey(rowFactor))
				clusterMap.put(rowFactor, new ArrayList<>());
			clusterMap.get(rowFactor).add(id);
		}
		return clusterMap;
	}

	public String getSortedRow() { return sortedRow; }

	public SortableTableModel getTableModel() {
		if (tableModel == null)
			tableModel = new FileCategoryTableModel(this);
		return tableModel;
	}

	public class FileCategoryTableModel extends SortableTableModel {
		final FileCategory design;
		final Experiment experiment;

		FileCategoryTableModel(final FileCategory design) {
			super(design.getHeaderCols());
			this.design = design;
			this.experiment = design.experiment;
			// NOTE: we're pivoting the table!
			// ncols = design.rows.size();
			// nrows = design.columns.length-1;
			hdrCols = 1;
		}

		@Override
		public int getColumnCount() { return design.getNCols(); }

		@Override
		public String getColumnName(int column) {
			if (columnIndex == null) 
				return strip(design.getColumnLabel(column));
			else
				return strip(design.getColumnLabel(columnIndex[column]));
		}

		@Override
		public int getRowCount() { 
			return design.getNRows();
		}

		@Override
		public Class<?> getColumnClass(int column) {
			return String.class;
		}

		@Override
		public Object getValueAt(int row, int column) {
			if (column == 0) {
				return strip(design.getRowLabel(row+1));
			}

			if (columnIndex != null)
				column = columnIndex[column];

			if (design.categories[row][column] == null)
				return "";

			return strip(design.categories[row][column]);
		}

		@Override
		public void sortColumns(int row) {
			sortedRow = strip(design.getRowLabel(row));
			super.sortColumns(row);
		}

		public String strip(String str) {
			return str.replaceAll("^\"|\"$", "");
		}
	}
}
