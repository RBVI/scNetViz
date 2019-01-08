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

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.AbstractCategory;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.SimpleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class GXACluster extends AbstractCategory implements IntegerMatrix {
	public static String GXA_CLUSTER_URI = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download?fileType=cluster";

	final Logger logger;

	// The suggested k value
	int k = 0;

	// The number of K values provided
	int nK = 0;

	// The lowest K value
	int minK = 0;

	// The clusters
	int[][] clusters;

	// K-sort
	int sortedK = -1;

	// The Table model
	GXAClusterTableModel tableModel = null;

	public GXACluster(final ScNVManager scManager, final GXAExperiment experiment) {
		super(scManager, experiment, "Cluster", 0, 0);
		super.hdrCols = 2;

		logger = Logger.getLogger(CyUserLog.NAME);
	}

	@Override
	public int getDefaultRow() { 
		if (k >= minK)
			return k-minK;
		return 0;
	}

	@Override
	public String toString() { return "Cluster";}

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
	public String mkLabel(Object cat) { return "Cluster "+cat.toString();}

	@Override
	public String getCategoryType() { return "Cluster";}

	@Override
	public Experiment getExperiment() { return experiment;}

	@Override
	public Matrix getMatrix() { return this; }

	@Override
	public String getMatrixType() { return "Simple Integer Matrix"; }

	@Override
	public int getIntegerValue(int row, int column) {
		return clusters[row][column];
	}

	@Override
	public int getIntegerValue(String row, String column) {
		int col = colLabels.indexOf(column);

		int intRow = Integer.valueOf(row);
		return clusters[intRow-minK][col-2];
	}

	@Override
	public int[][] getIntegerMatrix(int missing) { return clusters; }

	@Override
	public int getHeaderCols() { return 2; }

	public int getNRows() {
		return clusters[0].length-1;
	}

	public int[] getCluster() {
		if (k == 0)
			return getCluster(minK); // No default K was selected

		return getCluster(k);
	}

	public int[] getCluster(int kClust) {
		int ncolumns = clusters.length;
		int[] clustRow = new int[ncolumns];
		for (int col = 0; col < ncolumns; col++) {
			clustRow[col] = clusters[col][kClust-minK];
		}
		return clustRow;
	}

	public int getMinK() { return minK; }
	public int getK() { return k; }
	public int getSortedK() { return sortedK; }

	public Map<Integer,List<String>> getClusterList(int kClust) {
		Map<Integer, List<String>> clusterMap = new HashMap<>();
		int[] clusterArray = getCluster(kClust);
		for (int i = 0; i < clusterArray.length; i++) {
			int cluster = clusterArray[i];
			if (!clusterMap.containsKey(cluster)) {
				clusterMap.put(cluster, new ArrayList<>());
			}
			clusterMap.get(cluster).add(colLabels.get(i+2));
		}
		return clusterMap;
	}

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


	public static GXACluster fetchCluster(ScNVManager scManager, String accession, 
	                                      GXAExperiment experiment, TaskMonitor monitor) {
		// Get the URI
		List<String[]> input = CSVReader.readCSVFromHTTP(monitor, GXA_CLUSTER_URI, accession);
		if (input == null || input.size() < 2) return null;

		int nclusters = input.size();
		int ncolumns = input.get(0).length-2;

		GXACluster gxaCluster = new GXACluster(scManager, experiment);

		gxaCluster.clusters = new int[ncolumns][nclusters];
		boolean first = true;
		List<String> lbl = new ArrayList<>();

		for (String[] line: input) {
			if (first) {
				first = false;
				gxaCluster.setColLabels(Arrays.asList(line));
				// gxaCluster.headers = line;
				continue;
			}
			int thisK = Integer.parseInt(line[1]);
			lbl.add("k = "+thisK);
			if (line[0].equals("TRUE")) {
				gxaCluster.k = thisK;
			}
			if (gxaCluster.minK == 0)
				gxaCluster.minK = thisK;

			for (int i = 2; i < line.length; i++) {
				gxaCluster.clusters[i-2][thisK-gxaCluster.minK] = Integer.parseInt(line[i]);
			}

		}

		gxaCluster.setRowLabels(lbl);

		gxaCluster.selectedRow = gxaCluster.k - gxaCluster.minK;

		if (monitor != null) {
			monitor.showMessage(TaskMonitor.Level.INFO, "Read "+(nclusters-1)+" clusters.  Suggested K = "+gxaCluster.k);
			monitor.showMessage(TaskMonitor.Level.INFO, "   minimum K = "+gxaCluster.minK);
		} else {
			gxaCluster.logger.info("Read "+(nclusters-1)+" clusters.  Suggested K = "+gxaCluster.k);
			gxaCluster.logger.info("   minimum K = "+gxaCluster.minK);
		}

		gxaCluster.nCols = ncolumns;
		gxaCluster.nRows = nclusters;
		
		return gxaCluster;
	}

	public TableModel getTableModel() {
		if (tableModel == null)
			tableModel = new GXAClusterTableModel(this);
		return tableModel;
	}

	@Override
	public Object getValue(int row, int col) {
		return clusters[col][row];
	}

	public class GXAClusterTableModel extends SortableTableModel {
		final GXACluster cluster;
		final Experiment experiment;

		GXAClusterTableModel(final GXACluster cluster) {
			super(cluster.getHeaderCols());
			this.cluster = cluster;
			this.experiment = cluster.experiment;
			// System.out.println("ncols = "+ncols+", nrows = "+nrows);
			// System.out.println("clusters.length = "+cluster.clusters.length);
			hdrCols = 2;
		}

		@Override
		public int getColumnCount() { return cluster.getNCols()+2; }

		@Override
		public int getSelectedRow() { return cluster.getSelectedRow(); }

		@Override
		public void setSelectedRow(int selectedRow) { cluster.setSelectedRow(selectedRow); }

		@Override
		public String getColumnName(int column) {
			if (columnIndex != null)
				return cluster.getColumnLabel(columnIndex[column]);
			else
				return cluster.getColumnLabel(column);
		}

		@Override
		public int getRowCount() { 
			return cluster.getNRows();
		}

		@Override
		public Class getColumnClass(int column) {
			switch(column) {
				case 0:
					return String.class;
				default:
					return Integer.class;
			}
		}

		@Override
		public void sortColumns(int row) {
			sortedK = row+cluster.minK;
			super.sortColumns(row);
		}

		@Override
		public Object getValueAt(int row, int column) {
			// System.out.println("getValueAt: "+row+","+column);
			switch (column) {
				case 0:
					return (row+cluster.minK) == k ? "True" : "False";
				case 1:
					return row+cluster.minK;
				default:
					int value;
					if (columnIndex != null)
						value = cluster.clusters[columnIndex[column]-2][row];
					else
						value = cluster.clusters[column-2][row];
					return new Integer(value);
			}
		}
	}
}
