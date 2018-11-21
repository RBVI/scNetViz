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
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.StringMatrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.SimpleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.utils.LogUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class FileCategory extends SimpleMatrix implements Category {
	final Logger logger;
	final Experiment experiment;
	final String name;
	int hdrCols = 1;

	String[][] stringCategories = null;
	int[][] intCategories = null;
	double[][] doubleCategories = null;

	Map<Object, Integer> sizes = null;
	Map<Object, double[]> means = null;
	// double[][] means = null;

	final String dataType;

	String sortedRow = null;

	int selectedRow = -1;

	Map<Object, List<Integer>> catMap = null;

	SortableTableModel tableModel = null;

	public FileCategory(final ScNVManager scManager, 
	                    final Experiment experiment, final String name,
	                    final String type, int nRows, int nCols) {
		super(scManager);
		super.nRows = nRows;
		super.nCols = nCols;
		this.experiment = experiment;
		this.name = name;
		logger = Logger.getLogger(CyUserLog.NAME);
		this.dataType = type;
		if (dataType.equals("text"))
			stringCategories = new String[nCols][nRows];
		else if (dataType.equals("integer"))
			intCategories = new int[nCols][nRows];
		else if (dataType.equals("float"))
			doubleCategories = new double[nCols][nRows];
	}

	@Override
	public String toString() { return name;}

	@Override
	public String getCategoryType() { return name;}

	@Override
	public Experiment getExperiment() { return experiment;}

	@Override
	public String getMatrixType() { 
		return "Simple "+dataType+" matrix";
	}

	@Override
	public Matrix getMatrix() { return this;}

	@Override
	public int getHeaderCols() { return 1; }

	public void setValue(int row, int col, String value) {
		if (dataType.equals("text"))
			stringCategories[col][row] = value;
		else if (dataType.equals("integer"))
			intCategories[col][row] = Integer.parseInt(value);
		else if (dataType.equals("float"))
			doubleCategories[col][row] = Double.parseDouble(value);
	}

	public Object getValue(int row, int col) {
		if (dataType.equals("text"))
			return stringCategories[col][row];
		else if (dataType.equals("integer"))
			return new Integer(intCategories[col][row]);
		else if (dataType.equals("float"))
			return new Double(doubleCategories[col][row]);
		return null;
	}

	@Override
	public Map<Object, double[]> getMeans(int category) {
		if (means != null && category == selectedRow)
			return means;

		if (sizes == null || category != selectedRow) {
			getSizes(category);
		}

		means = new HashMap<>();

		Matrix mtx = experiment.getMatrix();
		DoubleMatrix dMat = null;
		IntegerMatrix iMat = null;
		if (mtx instanceof DoubleMatrix) {
			dMat = (DoubleMatrix)mtx;
		} else if (mtx instanceof IntegerMatrix) {
			iMat = (IntegerMatrix)mtx;
		}

		for (Object key: catMap.keySet()) {
			List<Integer> arrays = catMap.get(key);
			double[] catMean = new double[mtx.getNRows()];
			for (int row = 0; row < mtx.getNRows(); row++) {
				double mean = 0.0;
				for (Integer col: arrays) {
					if (dMat != null) {
						mean += dMat.getDoubleValue(row, col)/(double)arrays.size();
					} else if (iMat != null) {
						mean += 
							(double)iMat.getIntegerValue(row, col)/(double)arrays.size();
					}
				}
				catMean[row] = mean;
			}
			means.put(key, catMean);
		}
		return means;
	}

	@Override
	public Map<Object, Integer> getSizes(int category) {
		if (sizes == null || category != selectedRow) {
			// This creates the sizes map as a by-product
			getUniqValues(category);
		}
		return sizes;
	}

	// dDRthreshold is the cutoff for the minimum difference between clusters
	@Override
	public void filter(int category, double dDRthreshold) {
		return;
	}

	// Calculate the logGER between each category and all other categories
	// This will trigger the calculation of means and sizes
	@Override
	public Map<String, double[]> getLogGER(int category) {
		return null;
	}

	// Calculate the logGER between the category and all other categories
	// This will trigger the calculation of means and sizes
	@Override
	public double[] getLogGER(int category, String category1) {
		return null;
	};

	// Calculate the logGER between the two categories
	// This will trigger the calculation of means and sizes
	@Override
	public double[] getLogGER(int category, String category1, String category2) {
		return null;
	}

	public static FileCategory fetchCategory(ScNVManager scManager, Experiment experiment,
	                                         File file, String dataCategory, boolean transpose, int hdrCols,
	                                         TaskMonitor monitor) throws Exception {

		List<String[]> input = CSVReader.readCSV(monitor, file);
		if (input == null || input.size() < 2) {
			// System.out.println("No input!");
			return null;
		}

		int nRows = input.size()-1; // Rows don't include the header
		int nCols = input.get(0).length-hdrCols;
		if (transpose) {
			int x = nCols;
			nCols = nRows;
			nRows = x;
		}

		FileCategory fileCategory = new FileCategory(scManager, experiment, file.getName(), dataCategory, nRows, nCols);
		fileCategory.hdrCols = hdrCols;

		List<String> labels;

		if (!transpose) {
			fileCategory.setColLabels(Arrays.asList(input.get(0)));
			labels = new ArrayList<String>(fileCategory.nRows);
		} else {
			String[] colLabels = input.get(0);
			String[] newLabels = Arrays.copyOfRange(colLabels, 1, colLabels.length);
			fileCategory.setRowLabels(Arrays.asList(newLabels));
			labels = new ArrayList<String>(fileCategory.nCols);
			labels.add("Category");
		}


		boolean first = true;
		int lineNumber = 0;
		for (String[] line: input) {
			if (first) {
				first = false;
			} else {
				labels.add(line[0]);
				// System.out.println("Label["+(lineNumber-1)+"]: "+line[0]);
				if (!transpose) {
					for (int col = 0; col < fileCategory.nCols; col++) {
						fileCategory.setValue(lineNumber-1, col, line[col+1]);
					}
				} else {
					for (int row = 0; row < fileCategory.nRows; row++) {
						fileCategory.setValue(row, lineNumber-1, line[row+1]);
					}
				}
			}
			lineNumber++;
		}

		if (!transpose) {
			// System.out.println("Found "+labels.size()+" row labels");
			fileCategory.setRowLabels(labels);
		} else {
			// System.out.println("Found "+labels.size()+" column labels");
			fileCategory.setColLabels(labels);
		}

		LogUtils.log(monitor, TaskMonitor.Level.INFO, "Read "+fileCategory.nRows+
			                    " rows with "+fileCategory.nCols+" columns");
		// System.out.println("Read "+fileCategory.nRows+" rows with "+fileCategory.nCols+" columns");

		return fileCategory;
	}

	public String getSortedRow() { return sortedRow; }

	public int getSelectedRow() { return selectedRow; }
	public void setSelectedRow(int selectedRow) { this.selectedRow = selectedRow; }

	public SortableTableModel getTableModel() {
		if (tableModel == null)
			tableModel = new FileCategoryTableModel(this);
		return tableModel;
	}

	private int getUniqValues(int row) {
		catMap = new HashMap<>();
		sizes = new HashMap<>();
		for (int col = hdrCols; col < nCols; col++) {
			Object v = getValue(row, col);
			if (!catMap.containsKey(v)) {
				catMap.put(v, new ArrayList<>());
				sizes.put(v, -1);
			}
			catMap.get(v).add(col);
			sizes.put(v, sizes.get(v)+1);
		}
		return catMap.keySet().size();
	}

	public class FileCategoryTableModel extends SortableTableModel {
		final FileCategory category;
		final Experiment experiment;

		FileCategoryTableModel(final FileCategory category) {
			super(category.getHeaderCols());
			this.category = category;
			this.experiment = category.experiment;
			hdrCols = 1;
		}

		@Override
		public int getColumnCount() { return category.getNCols(); }

		@Override
		public String getColumnName(int column) {
			if (columnIndex == null) 
				return strip(category.getColumnLabel(column));
			else
				return strip(category.getColumnLabel(columnIndex[column]));
		}

		@Override
		public int getRowCount() { 
			return category.getNRows();
		}

		@Override
		public Class<?> getColumnClass(int column) {
			if (column < hdrCols)
				return String.class;
			if (category.dataType.equals("text"))
				return String.class;
			else if (category.dataType.equals("integer"))
				return Integer.class;
			else if (category.dataType.equals("double"))
				return Double.class;
			return String.class;
		}

		@Override
		public Object getValueAt(int row, int column) {
			if (column == 0) {
				return strip(category.getRowLabel(row));
			}

			if (columnIndex != null)
				column = columnIndex[column];

			if (category.getValue(row, column) == null)
				return "";

			if (category.dataType.equals("text"))
				return strip(category.getValue(row, column).toString());
			else 
				return category.getValue(row, column);
		}

		@Override
		public void sortColumns(int row) {
			sortedRow = strip(category.getRowLabel(row));
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
