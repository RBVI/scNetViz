package edu.ucsf.rbvi.scNetViz.internal.sources.hca;

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

import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
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

public class HCADesign extends AbstractCategory implements StringMatrix {
	public static String HCA_DESIGN_URI = "https://www.ebi.ac.uk/hca/sc/experiment/%s/download?fileType=experiment-design";
	final Logger logger;

	String[][] categories;

	String sortedRow = null;

	SortableTableModel tableModel = null;

	Source source;

	public HCADesign(final ScNVManager scManager, final HCAExperiment experiment) {
		super(scManager, experiment, "Design/Factors", 0, 0);
		logger = Logger.getLogger(CyUserLog.NAME);
		source = scManager.getSource("HCA");
		hdrCols = 1;
	}

	@Override
	public int getDefaultRow() { return -1;} // We don't have a default

	@Override
	public String toString() { return "Design/Factors";}

	@Override
	public String toJSON() { 
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		builder.append("\"name\": \""+toString()+"\",");
		builder.append("\"source\": \""+source+"\",");
		builder.append("\"rows\": "+getMatrix().getNRows()+",");
		builder.append("\"columns\": "+getMatrix().getNCols()+",");
		builder.append("\"default row\": "+getDefaultRow());
		builder.append("}");
		return builder.toString();
	}

	@Override
	public String getCategoryType() { return "Design/Factors";}

	@Override
	public Experiment getExperiment() { return experiment;}

	@Override
	public Source getSource() { return source; }

	@Override
	public String getMatrixType() { return "Simple String Matrix";}

	@Override
	public Matrix getMatrix() { return this;}

	@Override
	public String[][] getStringMatrix(boolean transpose, boolean excludeControls) { 
		return getStringMatrix(transpose);
	}

	@Override
	public String[][] getStringMatrix(boolean transpose) { 
		if (transpose && transposed)
			transpose = false;
		else
			transpose = transpose | transposed;
		if (!transpose)
			return categories;

		String[][] newArray = getStringMatrix(transpose);
		for (int row = 0; row < nRows; row++) {
			for (int col = 0; col < nCols; col++) {
				if (transposed)
					newArray[col][row] = categories[col][row];
				else
					newArray[row][col] = categories[row][col];
			}
		}
		return newArray;
	}

	@Override
	public String[][] getStringMatrix() { return categories;}

	@Override
	public String getValue(int row, int col) { 
		// System.out.println("categories["+row+"]["+col+"] = "+categories[row][col]);
		return categories[row][col];
	}

	@Override
	public String getValue(String rowLabel, String colLabel) { 
		int row = rowLabels.indexOf(rowLabel);
		int col = colLabels.indexOf(colLabel);
		return categories[row][col];
	}

	@Override
	public int getHeaderCols() { return hdrCols; }

	// dDRthreshold is the cutoff for the minimum difference between clusters
	@Override
	public void filter(int category, double dDRthreshold) {
		return;
	}

	public void fetchDesign(HCAExperiment experiment, List<String[]> inputTable, TaskMonitor monitor) {

		List<String[]> input = new ArrayList(inputTable);
		// For HCA, the design information is encoded in the cell table.
		if (input == null || input.size() < 2) return;

		nCols = input.size();
		nRows = input.get(0).length-1;
		// System.out.println("HCA fetchDesign: nCols = "+nCols+", nRows = "+nRows);

		setRowLabels(stripArray(input.get(0), 1));
		categories = new String[nRows][nCols];
		List<String> colLabels = new ArrayList<String>(nCols);
		colLabels.add("Category");

		boolean first = true;
		int col = hdrCols;
		for (String[] line: input) {
			if (col >= hdrCols) {
				if (first) {
					first = false;
				} else {
					String label = stripQuotes(new String(line[0]));
					colLabels.add(label);
				}
				for (int row = 1; row < nRows; row++) {
					categories[row-1][col-hdrCols] = stripQuotes(new String(line[row]));
				}
			}
			col++;

		}
		// System.out.println("colLabels.size() = "+colLabels.size());
		setColLabels(colLabels);

		LogUtils.log(monitor, TaskMonitor.Level.INFO, "Read "+nRows+
			                    " rows with "+nCols+" columns");
		return;
	}

	static private List<String> stripArray(String[] array, int offset) {
		List<String> result = new ArrayList<>();
		for (int i = offset; i < array.length; i++) {
			String str = new String(array[i]);
			result.add(str.replaceAll("^\"|\"$", ""));
		}
		return result;
	}

	static private String stripQuotes(String str) {
		return str.replaceAll("^\"|\"$", "");
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
			tableModel = new HCADesignTableModel(this);
		return tableModel;
	}

	public class HCADesignTableModel extends SortableTableModel {
		final HCADesign design;
		final Experiment experiment;

		HCADesignTableModel(final HCADesign design) {
			super(design.getHeaderCols());
			this.design = design;
			this.experiment = design.experiment;
			// NOTE: we're pivoting the table!
			// ncols = design.rows.size();
			// nrows = design.columns.length-1;
			hdrCols = 0;
		}

		@Override
		public int getColumnCount() { return design.getNCols(); }

		@Override
		public int getSelectedRow() { return design.getSelectedRow(); }

		@Override
		public void setSelectedRow(int selectedRow) { design.setSelectedRow(selectedRow); }

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
				return strip(design.getRowLabel(row));
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
