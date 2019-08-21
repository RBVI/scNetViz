package edu.ucsf.rbvi.scNetViz.internal.sources.hca;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.table.TableModel;

import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;

import org.apache.log4j.Logger;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVWriter;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;


public class HCAExperiment implements Experiment {
	public static String HCA_MATRIX_URL = "https://matrix.data.humancellatlas.org/v1/matrix/";

	final Logger logger;

	String accession = null;
	List<String[]> rowTable = null;
	List<String[]> colTable = null;
	MatrixMarket mtx = null;
	HCAMetadata hcaMetadata = null;
	List<Category> categories;
	double[][] tSNE;
	String plotType = "t-SNE";

	final ScNVManager scNVManager;
	final HCAExperiment hcaExperiment;
	final HCASource source;
	HCAExperimentTableModel tableModel = null;

	DifferentialExpression diffExp = null;

	public HCAExperiment (ScNVManager manager, HCASource source, HCAMetadata entry) {
		this.scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);
		this.hcaExperiment = this;
		categories = new ArrayList<Category>();
		categories.add(new HCADesign(manager, this));
		this.source = source;
		this.hcaMetadata = entry;
		this.accession = (String)hcaMetadata.get(Metadata.ACCESSION);
	}

	public Matrix getMatrix() { return mtx; }
	public String getAccession() { return accession; }
	public String getSpecies() { return (String)hcaMetadata.get(Metadata.SPECIES); }

	public List<String[]> getColumnLabels() { return colTable; }
	public List<String[]> getRowLabels() { return rowTable; }

	// public HCACluster getClusters() { return hcaCluster; }
	// public HCADesign getDesign() { return hcaDesign; }

	public List<Category> getCategories() { return categories; }

	public Category getCategory(String categoryName) { 
		for (Category cat: categories) {
			if (cat.toString().equals(categoryName))
				return cat;
		}
		return null;
	}

	public void addCategory(Category c) { categories.add(c); }

	public Metadata getMetadata() { return hcaMetadata; }

	public Category getDefaultCategory() { return categories.get(0); }

	public Source getSource() { return source; }

	@Override
	public void setTSNE(double[][] tsne) {
		tSNE = tsne;
	}

	@Override
	public double[][] getTSNE() {
		return tSNE;
	}

	@Override
	public void setPlotType(String type) { this.plotType = type; }

	@Override
	public String getPlotType() { return plotType; }

	public DifferentialExpression getDiffExp() { return diffExp; }
	public void setDiffExp(DifferentialExpression de) { diffExp = de; }

	public TableModel getTableModel() { 
		if (tableModel == null)
			tableModel = new HCAExperimentTableModel(scNVManager, this);
		return tableModel;
	}

	public void fetchMTX (final TaskMonitor monitor) {
		// Create the JSON argument
		String query = "{\"filter\": {\"op\":\"=\",\"value\":\""+
		                 accession+"\",\"field\":\"project.provenance.document_id\"},\"format\":\"mtx\"}";
		// Get the URI
		try {
			CloseableHttpClient httpclient = HttpClients.createDefault();
			JSONObject jsonResponse = HTTPUtils.postJSON(HCA_MATRIX_URL, httpclient, query, monitor);

			String requestId = (String)jsonResponse.get("request_id");

			// Now that we have a request ID, interate through until the matrix server is ready
			int retries = 120;
			JSONObject status = null;
			while (retries > 0) {
				status = HTTPUtils.fetchJSON(HCA_MATRIX_URL+requestId, httpclient, monitor);
				if (status != null && status.containsKey("status") && 
				    !((String)status.get("status")).equals("In Progress")) {
					break;
				}
				Thread.sleep(5000); // 5 seconds
				retries--;
			}
			if (retries == 0 || status == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "Timeout waiting for matrix service");
				httpclient.close();
				return;
			}

			String statusValue = (String)status.get("status");
			if (!statusValue.equals("Complete")) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "Error getting matrix: "+statusValue);
				httpclient.close();
				return;
			}
			httpclient.close();

			CloseableHttpClient httpZipClient = HttpClients.createDefault();
			String matrixURL = (String)status.get("matrix_url");
			ZipInputStream zipStream = HTTPUtils.getZipStream(matrixURL, httpZipClient, monitor);
			ZipEntry entry;
			while ((entry = zipStream.getNextEntry()) != null) {
				String name = entry.getName();
				if (name.endsWith("cells.tsv.gz")) {
					colTable = CSVReader.readCSV(monitor, zipStream, name, 0);

					// First, create our design category
					((HCADesign)categories.get(0)).fetchDesign(this, colTable, monitor);

					// Now, we need to skip over the header
					colTable.remove(0);

					if (mtx != null) 
						mtx.setColumnTable(colTable, 1);
				} else if (name.endsWith("genes.tsv.gz")) {
					rowTable = CSVReader.readCSV(monitor, zipStream, name, 1);
					if (mtx != null) 
						mtx.setRowTable(rowTable, 0);
				} else if (name.endsWith(".mtx.gz")) {
					mtx = new MatrixMarket(scNVManager, null, null);
					if (rowTable != null)
						mtx.setRowTable(rowTable, 0);

					if (colTable != null)
						mtx.setColumnTable(colTable, 1);

					mtx.readMTX(monitor, zipStream, name);
				}
			}
			zipStream.close();
			httpZipClient.close();

		} catch (Exception e) {}
		scNVManager.addExperiment(accession, this);
	}

	public void fetchDesign (final TaskMonitor monitor) {
	}

	public String toHTML() {
		return hcaMetadata.toHTML();
	}

	public String toString() {
		return getAccession();
	}

	public String toJSON() {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		builder.append("\"source\": \""+getSource().toString()+"\",\n");
		builder.append("\"source name\": \""+getSource().getName()+"\",\n");
		builder.append("\"metadata\": "+hcaMetadata.toJSON()+",\n");
		builder.append("\"rows\": "+getMatrix().getNRows()+",");
		builder.append("\"columns\": "+getMatrix().getNCols()+",");
		List<Category> categories = getCategories();
		builder.append("\"categories\": [");
		for (Category cat: categories) {
			builder.append(cat.toJSON()+",\n");
		}
		return builder.substring(0, builder.length()-2)+"]}";
	}

	public void createSessionFiles(String accession, List<File> files) throws Exception {
		String tmpDir = System.getProperty("java.io.tmpdir");
		String expPrefix = source.getName()+"."+accession;
		try {
			// Save the Experiment file as an MTX
			File mtxFile = new File(tmpDir, URLEncoder.encode(expPrefix)+".mtx");
			mtx.saveFile(mtxFile);
			files.add(mtxFile);
			File mtxRowFile = new File(tmpDir, URLEncoder.encode(expPrefix)+".mtx_rows");
			CSVWriter.writeCSV(mtxRowFile, rowTable);
			files.add(mtxRowFile);

			File mtxColFile = new File(tmpDir, URLEncoder.encode(expPrefix)+".mtx_cols");
			CSVWriter.writeCSV(mtxColFile, colTable);
			files.add(mtxColFile);
		} catch (Exception e) {
			logger.error("Unable to save MTX data for "+accession+" in session: "+e.toString());
				e.printStackTrace();
			return;
		}

		// Save each Category as a CSV
		for (Category cat: categories) {
			try {
				String catPrefix = URLEncoder.encode(expPrefix+"."+cat.getSource().getName()+"."+cat.toString());
				File catFile = new File(tmpDir, catPrefix+".csv");
				cat.saveFile(catFile);
				files.add(catFile);
			} catch (Exception e) {
				logger.error("Unable to save categtory data for "+accession+" "+cat+" in session: "+e.toString());
				e.printStackTrace();
			}
		}

		// Save the current DiffExp as a CSV
		if (diffExp != null) {
			try {
				String dePrefix = URLEncoder.encode(expPrefix+".diffExp");
				File deFile = new File(tmpDir, dePrefix+".csv");
				diffExp.saveFile(deFile);
				files.add(deFile);
			} catch (Exception e) {
				logger.error("Unable to save differential expression results for "+diffExp.toString()+" in session: "+e.toString());
				e.printStackTrace();
			}
		}
	}

	/*
	class FetchClusterThread implements Runnable {
		@Override
		public void run() {
			categories.set(0, HCACluster.fetchCluster(scNVManager, accession, hcaExperiment, null));
		}
	}

	class FetchDesignThread implements Runnable {
		@Override
		public void run() {
			categories.set(1, HCADesign.fetchDesign(scNVManager, accession, hcaExperiment, null));
		}
	}
	*/
}
