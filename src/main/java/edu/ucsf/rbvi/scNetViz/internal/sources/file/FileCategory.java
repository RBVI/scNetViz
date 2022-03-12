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

  /**
   * NOTE: the number of rows and columns should always reflect the
   * underlying matrix -- that is, should not include column or row
   * headers.
   */
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
		builder.append("\"type\": \""+dataType+"\",");
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
	public int getHeaderCols() { return hdrCols; }

	public void setValue(int row, int col, String value) {
    // System.out.println("Setting value["+row+"]["+col+"] to "+value);
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
	                                         File file, String dataCategory, String name, 
                                           boolean transpose, int hdrCols,
	                                         int keyCol, boolean zeroRelative,
	                                         TaskMonitor monitor) throws Exception {

    try {
		  List<String[]> input = CSVReader.readCSV(monitor, file);
  		if (input == null || input.size() < 2) {
  			System.out.println("No input!");
  			return null;
  		}

      if (name == null) 
        name = file.getName();
  		return createCategory(scManager, experiment, name, dataCategory, input,
  		                      transpose, hdrCols, keyCol, zeroRelative, null, monitor);
    } catch (FileNotFoundException e) {
		  LogUtils.log(monitor, TaskMonitor.Level.ERROR, "File not found: "+file.getName());
    } catch (IOException e) {
		  LogUtils.log(monitor, TaskMonitor.Level.ERROR, "Unable to read file: "+
                   file.getName()+" ["+e.getMessage()+"]");
    }
    return null;
	}

  // External interface to create a category from a file.  Note that there is a bit of trickery here.  For
  // some uses (e.g. our remote cluster routines) we provide a simple list of label,values.  This makes
  // managing the tables and labels a bit tough.  What we do is pass a static label (rowLabels) that provides
  // the row labels (there is only one of them) and set hdrCols to 0.  The tricky part is that *if* we have
  // a hdrCol of 0 and we have row labels, then by definition we actually have at least 1 header column.  But,
  // the data isn't structured that we, so we tell the matrix that there is a header column, but we process
  // the data as if there isn't.
	public static FileCategory createCategory(ScNVManager scManager, Experiment experiment,
	                                          String name, String dataCategory, List<String[]> lines,
	                                          boolean transpose, int hdrCols, int keyCol,
	                                          boolean zeroRelative, String[] rowLabels, TaskMonitor monitor) {

		int nRows = lines.size()-1; // Rows don't include the header
		int nCols = lines.get(0).length;
    System.out.println("nCols = "+nCols+", nRows = "+nRows);

		if (transpose) {
			int x = nCols;
			// Is this the matrix nCols?  If so, it should just be nRows
      // nCols = nRows+hdrCols;
      nCols = nRows;
      if (hdrCols == 0 && rowLabels != null)
        hdrCols = 1;
			nRows = x-1;
		} else {
      // Fixup the hdrCols
      if (hdrCols == 0 && rowLabels != null) {
        hdrCols = 1;
        nCols = nCols-1;
      } else if (rowLabels != null) {
        nCols = nCols-hdrCols;
      } else {
        if (hdrCols == 0)
          hdrCols = 1;
        nCols = nCols-hdrCols;
      }
    }
    System.out.println("Adjusted: nCols = "+nCols+", nRows = "+nRows);
    System.out.println("Adjusted: hdrCols = "+hdrCols+", rowLabels = "+rowLabels);

		FileCategory fileCategory = new FileCategory(scManager, experiment, name, dataCategory, nRows, nCols);
    fileCategory.setHdrCols(hdrCols);
		fileCategory.setHdrRows(1);  // XXX: Do we need multiple row headers?
    fileCategory.setColKey(keyCol);
		fileCategory.zeroRelative = zeroRelative;

		if (!transpose) {
      System.out.println("!tranpose");
      for (int row = 0; row < fileCategory.getHdrRows(); row++) {
        if (keyCol == 0) {
          // System.out.println("colLabels.size = "+lines.get(row).length);
			    fileCategory.setColLabels(lines.get(row), row);
        } else {
          // Set the key column as the first column
          String[] lbls = lines.get(row);
          setColLabels(fileCategory, lbls, keyCol, row);
        }
      }
		} else {
      // Note that we're using the matrix view of hdrCols
      for (int col = 0; col < fileCategory.getHdrCols(); col++) {
        if (rowLabels != null) {
			    fileCategory.setRowLabels(rowLabels, col);
        } else {
          String[] colLabels = lines.get(col);
          String[] newLabels = Arrays.copyOfRange(colLabels, fileCategory.getHdrRows(), colLabels.length);
          for (int lbl = 0; lbl < newLabels.length; lbl++)
            newLabels[lbl] = strip(newLabels[lbl]);
          fileCategory.setRowLabels(newLabels, col);
        }
      }
      fileCategory.setColLabel("Category", 0, 0);
		}

		boolean first = true;
		int lineNumber = 0;
		for (String[] line: lines) {
			if (first) {
        // Skip over the header row
				first = false;
			} else {
				if (!transpose) {
          setRowLabels(fileCategory, line, keyCol, lineNumber);
					for (int col = 0; col < fileCategory.nCols; col++) {
            try {
              fileCategory.setValue(lineNumber-1, col, line[col+hdrCols]);
            } catch (Exception e) {
		          LogUtils.log(monitor, TaskMonitor.Level.ERROR, 
                           "Unable to read value '"+line[col+hdrCols]+"' "+
                           "at "+(lineNumber-1)+","+col+" as type "+dataCategory);
              e.printStackTrace();
              return null;
            }
					}
				} else {
          if (keyCol == 0) {
            try {
            // System.out.println("Column["+lineNumber+"]="+strip(line[0]));
            fileCategory.setColLabel(strip(line), lineNumber);
            } catch (Exception e) { e.printStackTrace(); }
          } else {
            setColLabels(fileCategory, line, keyCol, lineNumber);
          }
					for (int row = 0; row < fileCategory.nRows; row++) {
            try {
              fileCategory.setValue(row, lineNumber-1, line[row+fileCategory.getHdrRows()]);
            } catch (Exception e) {
		          LogUtils.log(monitor, TaskMonitor.Level.ERROR, 
                           "Unable to read value '"+line[row+fileCategory.getHdrRows()]+"' "+
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

  private static void setColLabels(FileCategory fileCategory, String[] labels, int keyCol, int row) {
    // Set the key column as the first column
    int column = 1;
    fileCategory.setColLabel(strip(labels[keyCol]), row, 0);
    for (int col = 0; col < labels.length; col++) {
      if (col == keyCol) continue;
      fileCategory.setColLabel(strip(labels[col]), row, column);
      column++;
    }
  }

  private static void setRowLabels(FileCategory fileCategory, String[] labels, int keyCol, int row) {
    // Set the key column as the first column
    int column = 1;
    row = row-fileCategory.getHdrRows();
    fileCategory.setRowLabel(strip(labels[keyCol]), row, 0);
    for (int col = 0; col < fileCategory.getHdrCols(); col++) {
      if (col == keyCol) continue;
      fileCategory.setRowLabel(strip(labels[col]), row, column);
      column++;
    }
  }

  public static String strip(String str) {
    return str.replaceAll("^\"|\"$", "");
  }

  public static String[] strip(String[] str) {
    for (int i = 0; i < str.length; i++) {
      str[i] =  str[i].replaceAll("^\"|\"$", "");
    }
    return str;
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
			super(category.getHdrCols());
			this.category = category;
			this.experiment = category.experiment;
		}

		@Override
		public int getColumnCount() { return category.getNCols()+category.getHdrCols(); }

		@Override
		public String getColumnName(int column) {
			if (columnIndex == null) 
				return FileCategory.strip(category.getColumnLabel(column));
			else
				return FileCategory.strip(category.getColumnLabel(columnIndex[column]));
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
				return FileCategory.strip(category.getRowLabel(row, column));
			}

			if (columnIndex != null)
				column = columnIndex[column];

			if (category.getValue(row, column-hdrCols) == null) {
				return "";
      }

			if (category.dataType.equals("text") || dataType.equals("string")) {
				return FileCategory.strip(category.getValue(row, column-hdrCols).toString());
      } else  {
				return category.getValue(row, column-hdrCols);
      }
		}

		@Override
		public void sortColumns(int row) {
			sortedRow = FileCategory.strip(category.getRowLabel(row));
			super.sortColumns(row);
		}

		@Override
		public int getSelectedRow() { return category.getSelectedRow(); }

		@Override
		public void setSelectedRow(int selectedRow) { category.setSelectedRow(selectedRow); }
	}
}
