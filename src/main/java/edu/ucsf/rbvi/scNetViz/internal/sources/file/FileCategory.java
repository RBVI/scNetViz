package edu.ucsf.rbvi.scNetViz.internal.sources.file;

import java.io.File;
import java.io.InputStream;
import java.io.IOException;
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

import edu.ucsf.rbvi.scNetViz.internal.api.AbstractCategory;
import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.api.StringMatrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.SimpleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.utils.LogUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class FileCategory extends AbstractCategory implements Category {
	final Logger logger;

	String[][] stringCategories = null;
	int[][] intCategories = null;
	double[][] doubleCategories = null;

	final String dataType;

	String sortedRow = null;

	SortableTableModel tableModel = null;

	Source source = null;

	public FileCategory(final ScNVManager scManager, 
	                    final Experiment experiment, final String name,
	                    final String type, int nRows, int nCols) {
		super(scManager, experiment, name, nRows, nCols);
		logger = Logger.getLogger(CyUserLog.NAME);
		this.dataType = type;
		if (dataType.equals("text"))
			stringCategories = new String[nCols][nRows];
		else if (dataType.equals("integer"))
			intCategories = new int[nCols][nRows];
		else if (dataType.equals("float"))
			doubleCategories = new double[nCols][nRows];

		source = scManager.getSource("File");
	}

	@Override
	public String toString() { return name;}

	@Override
	public String toJSON() { 
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		builder.append("\"name\": \""+toString()+"\",");
		builder.append("\"source\": \""+source+"\",");
		builder.append("\"source name\": \""+source.getName()+"\",");
		builder.append("\"rows\": "+getMatrix().getNRows()+",");
		builder.append("\"columns\": "+getMatrix().getNCols()+",");
		builder.append("\"default row\": "+getDefaultRow());
		builder.append("}");
		return builder.toString();
	}

	@Override
	public String getCategoryType() { return name;}

	@Override
	public int getDefaultRow() { return -1;}

	@Override
	public Experiment getExperiment() { return experiment;}

	@Override
	public String getMatrixType() { 
		return "Simple "+dataType+" matrix";
	}

	@Override
	public Class<?> getMatrixClass() { 
		if (dataType.equals("float"))
			return Double.class;

		if (dataType.equals("integer"))
			return Integer.class;

		return String.class;
	}

	@Override
	public Matrix getMatrix() { return this;}

	@Override
	public Source getSource() { return source;}

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

	// dDRthreshold is the cutoff for the minimum difference between clusters
	@Override
	public void filter(int category, double dDRthreshold) {
		return;
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

	public SortableTableModel getTableModel() {
		if (tableModel == null)
			tableModel = new FileCategoryTableModel(this);
		return tableModel;
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
