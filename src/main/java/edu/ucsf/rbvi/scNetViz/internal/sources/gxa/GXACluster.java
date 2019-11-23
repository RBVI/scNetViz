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
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.SimpleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVWriter;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class GXACluster extends AbstractCategory implements IntegerMatrix {
	public static String GXA_CLUSTER_URI = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download?fileType=cluster";

	final Logger logger;

	// The suggested k value
	int suggestedK = 0;

	// The default row (corresponds to suggestedK)
	int defaultRow = 0;

	// The number of K values provided
	int nK = 0;

	// The clusters
	int[][] clusters;

	// K-sort
	int sortedK = -1;

	// The Table model
	GXAClusterTableModel tableModel = null;

	Source source;

	public GXACluster(final ScNVManager scManager, final GXAExperiment experiment, List<String[]> table) {
		this(scManager, experiment);

	}

	public GXACluster(final ScNVManager scManager, final GXAExperiment experiment) {
		super(scManager, experiment, "Cluster", 0, 0);
		super.hdrCols = 2;

		logger = Logger.getLogger(CyUserLog.NAME);
		source = scManager.getSource("GXA");
	}

	@Override
	public int getDefaultRow() { 
		return defaultRow;
	}

	@Override
	public String toString() { return "Cluster";}

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
		builder.append("\"suggested K\": "+suggestedK);
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
	public Source getSource() { return source;}

	@Override
	public Matrix getMatrix() { return this; }

	@Override
	public String getMatrixType() { return "Simple Integer Matrix"; }

	@Override
	public int getIntegerValue(int row, int column) {
		return clusters[column][row];
	}

	@Override
	public int getIntegerValue(String row, String column) {
		int col = colLabels.indexOf(column);

		int intRow = Integer.valueOf(row);
		return clusters[col-2][intRow];
	}

	@Override
	public int[][] getIntegerMatrix(int missing) { return clusters; }

	@Override
	public int[][] getIntegerMatrix(int missing, boolean transpose, boolean excludeControls) { 
		// We don't have controls, so just ignore it
		return getIntegerMatrix(missing, transpose);
	}

	@Override
	public int[][] getIntegerMatrix(int missing, boolean transpose) { 
		if (transpose && transposed)
			transpose = false;
		else
			transpose = transpose | transposed;
		if (!transpose)
			return clusters;

		int[][] newArray = getIntegerMatrix(transpose);
		for (int row = 0; row < nRows; row++) {
			for (int col = 0; col < nCols; col++) {
				if (transposed)
					newArray[col][row] = clusters[col][row];
				else
					newArray[row][col] = clusters[row][col];
			}
		}
		return newArray;
	}

	@Override
	public int getHeaderCols() { return hdrCols; }

	public int getNRows() {
		return clusters[0].length-1;
	}

	public int[] getCluster() {
		return getCluster(suggestedK);
	}

	public int[] getCluster(int kClust) {
		int ncolumns = clusters.length;
		int[] clustRow = new int[ncolumns];
		int kRow = getRowForK(kClust);
		for (int col = 0; col < ncolumns; col++) {
			clustRow[col] = clusters[col][kRow];
		}
		return clustRow;
	}

	public int getK() { return suggestedK; }
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

	public int getRowForK(int k) {
		for (int row = 0; row < nRows; row++) {
			if ((Integer)getValue(row, 0) == k)
				return row;
		}
		return -1;
	}

	public static GXACluster readCluster(ScNVManager scManager, GXAExperiment experiment, File file, 
	                                     JSONObject jsonCategory) throws IOException {
		List<String[]> input = CSVReader.readCSV(null, file);
		if (input == null || input.size() < 2) return null;

		GXACluster cluster =  getClusterFromCSV(scManager, experiment, input, null);
		if (jsonCategory.containsKey("suggested K")) {
			cluster.suggestedK = ((Long)jsonCategory.get("suggested K")).intValue();
			cluster.selectedRow = cluster.getRowForK(cluster.suggestedK);
		}
		return cluster;
	}

	public static GXACluster fetchCluster(ScNVManager scManager, String accession, 
	                                      GXAExperiment experiment, TaskMonitor monitor) {
		// Get the URI
		List<String[]> input = CSVReader.readCSVFromHTTP(monitor, GXA_CLUSTER_URI, accession);
		if (input == null || input.size() < 2) return null;

		return getClusterFromCSV(scManager, experiment, input, monitor);
	}

	private static GXACluster getClusterFromCSV(ScNVManager scManager, GXAExperiment experiment, 
	                                            List<String[]> input, TaskMonitor monitor) {
		int nclusters = input.size();
		int ncolumns = input.get(0).length;

		GXACluster gxaCluster = new GXACluster(scManager, experiment);

		// System.out.println("ncolumns = "+ncolumns+", ncluster = "+nclusters);
		gxaCluster.clusters = new int[ncolumns][nclusters];
		boolean first = true;
		List<String> lbl = new ArrayList<>();

		int clustering = 0;
		for (String[] line: input) {
			if (first) {
				first = false;
				// System.out.println("Column header line has: "+line.length+" columns");
				gxaCluster.setColLabels(Arrays.asList(line));
				// gxaCluster.headers = line;
				continue;
			}
			// System.out.println("Input:");
			// for (String l: line) { System.out.print(l+","); }
			// System.out.println();
			int thisK = Integer.parseInt(line[1]);
			lbl.add("k = "+thisK);
			if (line[0].equalsIgnoreCase("TRUE")) {
				gxaCluster.suggestedK = thisK;
				gxaCluster.selectedRow = clustering;
			}

			// System.out.println("line.length = "+line.length);
			// System.out.println("thisK = "+thisK);

			gxaCluster.clusters[0][clustering] = thisK;
			for (int i = 2; i < line.length; i++) {
				gxaCluster.clusters[i-1][clustering] = Integer.parseInt(line[i]);
			}
			clustering++;

		}

		gxaCluster.setRowLabels(lbl);

		if (monitor != null) {
			monitor.showMessage(TaskMonitor.Level.INFO, "Read "+(nclusters-1)+" clusters.  Suggested K = "+gxaCluster.suggestedK);
		} else {
			gxaCluster.logger.info("Read "+(nclusters-1)+" clusters.  Suggested K = "+gxaCluster.suggestedK);
		}

		gxaCluster.nCols = ncolumns;
		gxaCluster.nRows = nclusters;

		gxaCluster.source = experiment.getSource();
		
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
			// System.out.println("ncols = "+cluster.nCols+", nrows = "+cluster.nRows);
			// System.out.println("clusters.length = "+cluster.clusters.length);
			hdrCols = 2;
		}

		@Override
		public int getColumnCount() { return cluster.getNCols(); }

		@Override
		public int getSelectedRow() { return cluster.getSelectedRow(); }

		@Override
		public void setSelectedRow(int selectedRow) { 
			cluster.setSelectedRow(selectedRow); 
		}

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
			sortedK = row;
			super.sortColumns(row);
		}

		@Override
		public Object getValueAt(int row, int column) {
			// System.out.println("getValueAt: "+row+","+column);
			switch (column) {
				case 0:
					int thisK = (Integer) getValueAt(row, 1);
					return (thisK) == suggestedK ? "True" : "False";
				default:
					int value;
					if (columnIndex != null)
						value = cluster.clusters[columnIndex[column]-1][row];
					else
						value = cluster.clusters[column-1][row];
					return new Integer(value);
			}
		}
	}
}
