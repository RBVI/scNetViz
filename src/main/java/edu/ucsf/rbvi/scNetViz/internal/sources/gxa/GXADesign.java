package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.io.InputStream;
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

	public GXADesign(final ScNVManager scManager, final GXAExperiment experiment) {
		super(scManager, experiment, "Design/Factors", 0, 0);
		logger = Logger.getLogger(CyUserLog.NAME);
	}

	@Override
	public int getDefaultRow() { return -1;} // We don't have a default

	@Override
	public String toString() { return "Design/Factors";}

	@Override
	public String toJSON() { 
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		builder.append("name: "+toString()+",");
		builder.append("rows: "+getMatrix().getNRows()+",");
		builder.append("columns: "+getMatrix().getNCols()+",");
		builder.append("default row: "+getDefaultRow()+",");
		builder.append("}");
		return builder.toString();
	}

	@Override
	public String getCategoryType() { return "Design/Factors";}

	@Override
	public Experiment getExperiment() { return experiment;}

	@Override
	public String getMatrixType() { return "Simple String Matrix";}

	@Override
	public Matrix getMatrix() { return this;}

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
	public int getHeaderCols() { return 1; }

	// dDRthreshold is the cutoff for the minimum difference between clusters
	@Override
	public void filter(int category, double dDRthreshold) {
		return;
	}

	/*
	// Calculate the logGER between each category and all other categories
	// This will trigger the calculation of means and sizes
	@Override
	public Map<Object, double[]> getLogGER(int category) {
		return null;
	}

	// Calculate the logGER between the category and all other categories
	// This will trigger the calculation of means and sizes
	@Override
	public double[] getLogGER(int category, Object category1) {
		return null;
	};

	// Calculate the logGER between the two categories
	// This will trigger the calculation of means and sizes
	@Override
	public double[] getLogGER(int category, Object category1, Object category2) {
		return null;
	}
	*/

	public static GXADesign fetchDesign(ScNVManager scManager, String accession, 
	                                    GXAExperiment experiment, TaskMonitor monitor) {
		// Get the URI
		List<String[]> input = CSVReader.readCSVFromHTTP(monitor, GXA_DESIGN_URI, accession);
		if (input == null || input.size() < 2) return null;

		GXADesign gxaDesign = new GXADesign(scManager, experiment);
		gxaDesign.nCols = input.size()-1;
		gxaDesign.nRows = input.get(0).length-1;

		// System.out.println("nCols = "+gxaDesign.nCols+", experiment ncols = "+experiment.getMatrix().getNCols());

		gxaDesign.setRowLabels(stripArray(input.get(0), 1));
		gxaDesign.categories = new String[gxaDesign.nRows][gxaDesign.nCols];
		List<String> colLabels = new ArrayList<String>(gxaDesign.nCols);
		colLabels.add("Category");

		boolean first = true;
		int col = 0;
		for (String[] line: input) {
			if (first) {
				first = false;
			} else {
				colLabels.add(stripQuotes(line[0]));
				for (int row = 1; row < gxaDesign.nRows; row++) {
					gxaDesign.categories[row-1][col] = stripQuotes(line[row]);
				}
				col++;
			}
		}
		gxaDesign.setColLabels(colLabels);

		LogUtils.log(monitor, TaskMonitor.Level.INFO, "Read "+gxaDesign.nRows+
			                    " rows with "+gxaDesign.nCols+" columns");
		return gxaDesign;
	}

	static private List<String> stripArray(String[] array, int offset) {
		List<String> result = new ArrayList<>();
		for (int i = offset; i < array.length; i++) {
			String str = array[i];
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
