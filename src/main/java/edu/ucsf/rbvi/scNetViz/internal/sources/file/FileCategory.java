package edu.ucsf.rbvi.scNetViz.internal.sources.file;

import java.io.File;
import java.io.FileNotFoundException;
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

	boolean zeroRelative = false;

	public FileCategory(final ScNVManager scManager, 
	                    final Experiment experiment, final String name,
	                    final String type, int nRows, int nCols) {
		super(scManager, experiment, name, nRows, nCols);
		logger = Logger.getLogger(CyUserLog.NAME);
		this.dataType = type;
		if (dataType.equals("text") || dataType.equals("string"))
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
	public int getDefaultRow() { 
		if (nRows == 1) return 0;
		return -1;
	}

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
    try {
  		if (dataType.equals("text") || dataType.equals("string"))
  			stringCategories[col][row] = value;
  		else if (dataType.equals("integer")) {
  			intCategories[col][row] = Integer.parseInt(value);
  			if (zeroRelative)
  				intCategories[col][row]++;
  		} else if (dataType.equals("float"))
  			doubleCategories[col][row] = Double.parseDouble(value);
    } catch (Exception e) {
		  LogUtils.error("Unable to read value '"+value+"' as type "+dataType);
      return;
    }
	}

	public Object getValue(int row, int col) {
		if (dataType.equals("text") || dataType.equals("string"))
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
	                                         boolean zeroRelative,
	                                         TaskMonitor monitor) throws Exception {

    try {
		  List<String[]> input = CSVReader.readCSV(monitor, file);
  		if (input == null || input.size() < 2) {
  			// System.out.println("No input!");
  			return null;
  		}
  		return createCategory(scManager, experiment, file.getName(), dataCategory, input,
  		                      transpose, hdrCols, zeroRelative, monitor);
    } catch (FileNotFoundException e) {
		  LogUtils.log(monitor, TaskMonitor.Level.ERROR, "File not found: "+file.getName());
    } catch (IOException e) {
		  LogUtils.log(monitor, TaskMonitor.Level.ERROR, "Unable to read file: "+
                   file.getName()+" ["+e.getMessage()+"]");
    }
    return null;
	}

	public static FileCategory createCategory(ScNVManager scManager, Experiment experiment,
	                                          String name, String dataCategory, List<String[]> lines,
	                                          boolean transpose, int hdrCols,
	                                          boolean zeroRelative, TaskMonitor monitor) {

		int nRows = lines.size()-1; // Rows don't include the header
		int nCols = lines.get(0).length;
		if (transpose) {
			int x = nCols;
			nCols = nRows;
			nRows = x;
		}

		System.out.println("nCols = "+nCols+", nRows = "+nRows);

		FileCategory fileCategory = new FileCategory(scManager, experiment, name, dataCategory, nRows, nCols);
		fileCategory.setHdrCols(hdrCols);
		fileCategory.setHdrRows(1);  // XXX: Do we need multiple row headers?
		fileCategory.zeroRelative = zeroRelative;

		if (!transpose) {
      for (int row = 0; row < fileCategory.getHdrRows(); row++) {
			  fileCategory.setColLabels(lines.get(row), row);
      }
		} else {
      for (int col = 0; col < fileCategory.getHdrCols(); col++) {
			  String[] colLabels = lines.get(col);
			  String[] newLabels = Arrays.copyOfRange(colLabels, fileCategory.getHdrRows(), colLabels.length);
			  fileCategory.setRowLabels(newLabels, col);
      }
      fileCategory.setRowLabel("Category", 0, 0);
		}


		boolean first = true;
		int lineNumber = 0;
		for (String[] line: lines) {
			if (first) {
				first = false;
			} else {
				if (!transpose) {
          fileCategory.setRowLabel(line, lineNumber);
					for (int col = 0; col < fileCategory.nCols-hdrCols; col++) {
            try {
              System.out.println("Setting value for "+line[0]+","+fileCategory.getColumnLabel(col+hdrCols)+" to "+line[col+hdrCols]);
              fileCategory.setValue(lineNumber-1, col, line[col+hdrCols]);
            } catch (Exception e) {
		          LogUtils.log(monitor, TaskMonitor.Level.ERROR, 
                           "Unable to read value '"+line[col+hdrCols]+"' "+
                           "at "+lineNumber+","+col+" as type "+dataCategory);
              return null;
            }
					}
				} else {
          fileCategory.setColLabel(line, lineNumber);
					for (int row = 0; row < fileCategory.nRows; row++) {
            try {
						  fileCategory.setValue(row, lineNumber-1, line[row+hdrCols]);
            } catch (Exception e) {
		          LogUtils.log(monitor, TaskMonitor.Level.ERROR, 
                           "Unable to read value '"+line[row+hdrCols]+"' "+
                           "at "+lineNumber+","+row+" as type "+dataCategory);
              return null;
            }
					}
				}
			}
			lineNumber++;
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
			// hdrCols = 1;
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
			if (category.dataType.equals("text") || dataType.equals("string"))
				return String.class;
			else if (category.dataType.equals("integer"))
				return Integer.class;
			else if (category.dataType.equals("double"))
				return Double.class;
			return String.class;
		}

		@Override
		public Object getValueAt(int row, int column) {
			if (column < hdrCols) {
				return strip(category.getRowLabel(row));
			}

			if (columnIndex != null)
				column = columnIndex[column];

			if (category.getValue(row, column-hdrCols) == null)
				return "";

			if (category.dataType.equals("text") || dataType.equals("string"))
				return strip(category.getValue(row, column-hdrCols).toString());
			else 
				return category.getValue(row, column-hdrCols);
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
