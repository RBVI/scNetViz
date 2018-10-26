package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.io.InputStream;
import java.util.ArrayList;
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

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;

public class GXADesign implements Category {
	public static String GXA_DESIGN_URI = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download?fileType=experiment-design";
	final Logger logger;
	final ScNVManager scManager;
	final GXAExperiment experiment;

	// GXADesignTableModel tableModel = null;

	// The clusters
	String[] columns;
	List<String[]> rows;

	String sortedRow = null;

	public GXADesign(final ScNVManager scManager, final GXAExperiment experiment) {
		this.scManager = scManager;
		this.experiment = experiment;
		logger = Logger.getLogger(CyUserLog.NAME);
	}

	@Override
	public String getCategoryType() { return "Design/Factors";}

	public static GXADesign fetchDesign(ScNVManager scManager, String accession, 
	                                    GXAExperiment experiment, TaskMonitor monitor) {
		// Get the URI
		List<String[]> input = CSVReader.readCSVFromHTTP(monitor, GXA_DESIGN_URI, accession);
		if (input == null || input.size() < 2) return null;

		int nrows = input.size();
		int ncolumns = input.get(0).length;

		GXADesign gxaDesign = new GXADesign(scManager, experiment);

		gxaDesign.columns = input.get(0);
		gxaDesign.rows = new ArrayList<>();

		boolean first = true;
		for (String[] line: input) {
			if (first) {
				first = false;
			} else
				gxaDesign.rows.add(line);
		}

		if (monitor != null) {
			monitor.showMessage(TaskMonitor.Level.INFO, "Read "+gxaDesign.rows.size()+
			                    " rows with "+gxaDesign.columns.length+" columns");
		} else {
			gxaDesign.logger.info("Read "+gxaDesign.rows.size()+
			                      " rows with "+gxaDesign.columns.length+" columns");
		}

		return gxaDesign;
	}

	public Map<String,List<String>> getClusterList(String factor) {
		Map<String, List<String>> clusterMap = new HashMap<>();

		int factorColumn;
		for (factorColumn = 0; factorColumn < columns.length; factorColumn++) {
			if (factor.equals(columns[factorColumn]))
				break;
		}

		for (String[] row: rows) {
			String id = row[0];
			String rowFactor = row[factorColumn];
			if (!clusterMap.containsKey(rowFactor))
				clusterMap.put(rowFactor, new ArrayList<>());
			clusterMap.get(rowFactor).add(id);
		}
		return clusterMap;
	}

	public String getSortedRow() { return sortedRow; }

	/*
	public GXADesignTableModel getTableModel() {
		if (tableModel == null)
			tableModel = new GXADesignTableModel(this);
		return tableModel;
	}

	public class GXADesignTableModel extends GXASubTableModel {
		final GXADesign design;
		final GXAExperiment experiment;

		GXADesignTableModel(final GXADesign design) {
			super();
			this.design = design;
			this.experiment = design.experiment;
			// NOTE: we're pivoting the table!
			ncols = design.rows.size();
			nrows = design.columns.length-1;
			hdrCols = 1;
		}

		@Override
		public int getColumnCount() { return ncols; }

		@Override
		public String getColumnName(int column) {
			if (column == 0) return "Characteristic/Factor";
			if (columnIndex == null) 
				return strip(design.rows.get(column-1)[0]);
			else
				return strip(design.rows.get(columnIndex[column]-1)[0]);
		}

		@Override
		public int getRowCount() { 
			return nrows;
		}

		@Override
		public Class<?> getColumnClass(int column) {
			return String.class;
		}

		@Override
		public Object getValueAt(int row, int column) {
			if (column == 0) {
				return strip(columns[row+1]);
			}

			if (columnIndex != null)
				column = columnIndex[column-1];

			String[] line = rows.get(column);
			return strip(line[row+1]);
		}

		@Override
		public void sortColumns(int row) {
			sortedRow = strip(columns[row+1]);
			super.sortColumns(row);
		}

		public String strip(String str) {
			return str.replaceAll("^\"|\"$", "");
		}
	}
	*/
}
