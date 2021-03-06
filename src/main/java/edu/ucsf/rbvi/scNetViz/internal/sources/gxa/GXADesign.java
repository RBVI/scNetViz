package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

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

import org.json.simple.JSONObject;

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

public class GXADesign extends AbstractCategory implements StringMatrix {
	public static String GXA_DESIGN_URI = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download?fileType=experiment-design";
	final Logger logger;

	String[][] categories;

	String sortedRow = null;

	SortableTableModel tableModel = null;

	Source source;

	public GXADesign(final ScNVManager scManager, final GXAExperiment experiment) {
		super(scManager, experiment, "Design/Factors", 0, 0);
		logger = Logger.getLogger(CyUserLog.NAME);
		source = scManager.getSource("GXA");
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
		builder.append("\"source name\": \""+source.getName()+"\",");
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
		return categories[row][col+1];
	}

	@Override
	public String getValue(String rowLabel, String colLabel) { 
		int row = rowLabels.indexOf(rowLabel);
		int col = colLabels.indexOf(colLabel);
		return categories[row][col-1];
	}

	@Override
	public int getHeaderCols() { return 1; }

	// dDRthreshold is the cutoff for the minimum difference between clusters
	@Override
	public void filter(int category, double dDRthreshold) {
		return;
	}

	public static GXADesign readDesign(ScNVManager scManager, GXAExperiment experiment, File file, 
	                                   JSONObject jsonCategory) throws IOException {
		List<String[]> input = CSVReader.readCSV(null, file);
		if (input == null || input.size() < 2) return null;

		return getDesignFromSession(scManager, experiment, input);
	}

	public static GXADesign fetchDesign(ScNVManager scManager, String accession, 
	                                    GXAExperiment experiment, TaskMonitor monitor) {
		// Get the URI
		List<String[]> input = CSVReader.readCSVFromHTTP(monitor, GXA_DESIGN_URI, accession);
		if (input == null || input.size() < 2) return null;

		return getDesignFromCSV(scManager, experiment, input, monitor);
	}

	/**
	 * This is just like getDesignFromCSV except we don't transpose since it's already transposed in the
	 * session file
	 */
	private static GXADesign getDesignFromSession(ScNVManager scManager, GXAExperiment experiment, List<String[]> input) {
		GXADesign gxaDesign = new GXADesign(scManager, experiment);
		gxaDesign.nRows = input.size()-1;
		gxaDesign.nCols = input.get(0).length;

		System.out.println("nCols = "+gxaDesign.nCols+", nRows = "+gxaDesign.nRows);

		gxaDesign.setColLabels(stripArray(input.get(0), 0));
    gxaDesign.categories = new String[gxaDesign.nRows][gxaDesign.nCols];
    // List<String> rowLabels = new ArrayList<String>(gxaDesign.nRows);

		boolean first = true;
		int row = 0;
		for (String[] line: input) {
			if (first) {
				first = false;
				// System.out.println("Column header line has: "+line.length+" columns");
				// gxaCluster.headers = line;
				continue;
			}
			// rowLabels.add(stripQuotes(line[0]));
      gxaDesign.setRowLabel(stripQuotes(line[0]), row, 0);
			for (int col = 1; col < gxaDesign.nCols+1; col++) {
				gxaDesign.categories[row][col-1] = stripQuotes(line[col]);
			}
			row++;
		}
		// gxaDesign.setRowLabels(rowLabels);
		return gxaDesign;
	}

	private static GXADesign getDesignFromCSV(ScNVManager scManager, GXAExperiment experiment, List<String[]> input, TaskMonitor monitor) {
		GXADesign gxaDesign = new GXADesign(scManager, experiment);
		gxaDesign.nCols = input.size();
		gxaDesign.nRows = input.get(0).length-1;

		// System.out.println("nCols = "+gxaDesign.nCols+", experiment ncols = "+experiment.getMatrix().getNCols());
		// System.out.println("nCols = "+gxaDesign.nCols+", nRows = "+gxaDesign.nRows);

		gxaDesign.setRowLabels(stripArray(input.get(0), 1));
		gxaDesign.categories = new String[gxaDesign.nRows][gxaDesign.nCols];
		gxaDesign.setColLabel("Category", 0, 0);

		boolean first = true;
		int col = 1;
		for (String[] line: input) {
			if (first) {
				first = false;
			} else {
				gxaDesign.setColLabel(stripQuotes(line[0]), 0, col);
				for (int row = 1; row < gxaDesign.nRows; row++) {
					gxaDesign.categories[row-1][col] = stripQuotes(line[row]);
				}
				col++;
			}
		}
		gxaDesign.source = experiment.getSource();

		LogUtils.log(monitor, TaskMonitor.Level.INFO, "Read "+gxaDesign.nRows+
			                    " rows with "+gxaDesign.nCols+" columns");
		return gxaDesign;
	}

	static private List<String[]> stripArray(String[] array, int offset) {
		List<String[]> result = new ArrayList<>();
		for (int i = offset; i < array.length; i++) {
			String[] str = new String[1];
      str[0] = array[i].replaceAll("^\"|\"$", "");
			result.add(str);
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
			String[] id = rowLabels.get(row);
			String rowFactor = categories[row][factorColumn];
			if (!clusterMap.containsKey(rowFactor))
				clusterMap.put(rowFactor, new ArrayList<>());
			clusterMap.get(rowFactor).add(id[0]);
		}
		return clusterMap;
	}

	public String getSortedRow() { return sortedRow; }

	public SortableTableModel getTableModel() {
		if (tableModel == null)
			tableModel = new GXADesignTableModel(this);
		return tableModel;
	}

	public class GXADesignTableModel extends SortableTableModel {
		final GXADesign design;
		final Experiment experiment;

		GXADesignTableModel(final GXADesign design) {
			super(design.getHeaderCols());
			this.design = design;
			this.experiment = design.experiment;
			// NOTE: we're pivoting the table!
			// ncols = design.rows.size();
			// nrows = design.columns.length-1;
			hdrCols = 1;
      hdrRows = 1;
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
