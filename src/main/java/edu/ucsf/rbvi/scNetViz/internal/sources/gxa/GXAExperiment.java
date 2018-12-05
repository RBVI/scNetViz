package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.table.TableModel;

import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
import org.apache.log4j.Logger;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;

public class GXAExperiment implements Experiment {
	public static String RESULTS_URL = "https://www.ebi.ac.uk/gxa/sc/experiments/%s/Results";
	public static String GXA_MTX_URI = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download/zip?fileType=quantification-filtered";
	final Logger logger;

	String accession = null;
	List<String[]> rowTable = null;
	List<String[]> colTable = null;
	MatrixMarket mtx = null;
	GXAMetadata gxaMetadata = null;
	List<Category> categories;
	// GXACluster gxaCluster = null;
	// GXAIDF gxaIDF = null;
	// GXADesign gxaDesign = null;

	final ScNVManager scNVManager;
	final GXAExperiment gxaExperiment;
	final GXASource source;

	public GXAExperiment (ScNVManager manager, GXASource source, GXAMetadata entry) {
		this.scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);
		this.gxaExperiment = this;
		categories = new ArrayList<Category>();
		categories.add(new GXACluster(manager, this));
		categories.add(new GXADesign(manager, this));
		this.source = source;
		this.gxaMetadata = entry;
		this.accession = (String)gxaMetadata.get(Metadata.ACCESSION);
	}

	public Matrix getMatrix() { return mtx; }
	public String getAccession() { return accession; }
	public String getSpecies() { return (String)gxaMetadata.get(Metadata.SPECIES); }

	public List<String[]> getColumnLabels() { return colTable; }
	public List<String[]> getRowLabels() { return rowTable; }

	// public GXACluster getClusters() { return gxaCluster; }
	// public GXADesign getDesign() { return gxaDesign; }

	public List<Category> getCategories() { return categories; }

	public void addCategory(Category c) { categories.add(c); }

	public Metadata getMetadata() { return gxaMetadata; }

	public Category getDefaultCategory() { return categories.get(0); }

	public Source getSource() { return source; }

	public TableModel getTableModel() { return new GXAExperimentTableModel(scNVManager, this); }

	public void fetchMTX (final TaskMonitor monitor) {
		// Get the URI
		String fetchString = String.format(GXA_MTX_URI, accession);

		try {
			CloseableHttpClient httpclient = HttpClients.createDefault();
			HttpGet httpGet = new HttpGet(fetchString);
			CloseableHttpResponse response1 = httpclient.execute(httpGet);
			if (response1.getStatusLine().getStatusCode() != 200) {
				return;
			}
			HttpEntity entity1 = response1.getEntity();

			try {
				ZipInputStream zipStream = new ZipInputStream(entity1.getContent());

				ZipEntry entry;
				while ((entry = zipStream.getNextEntry()) != null) {
					String name = entry.getName();
					System.out.println("Name = "+name);
					if (name.endsWith(".mtx_cols")) {
						colTable = CSVReader.readCSV(monitor, zipStream, name);
						if (mtx != null) 
							mtx.setColumnTable(colTable);
					} else if (name.endsWith(".mtx_rows")) {
						rowTable = CSVReader.readCSV(monitor, zipStream, name);
						if (mtx != null) 
							mtx.setRowTable(rowTable);
					} else if (name.endsWith(".mtx")) {
						mtx = new MatrixMarket(scNVManager, null, null);
						mtx.setRowTable(rowTable);
						mtx.setColumnTable(colTable);
						mtx.readMTX(monitor, zipStream, name);
					}
					zipStream.closeEntry();
				}
				zipStream.close();
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				response1.close();
			}
		} catch (Exception e) {}
		scNVManager.addExperiment(accession, this);
	}

	public void fetchClusters (final TaskMonitor monitor) {
		categories.set(0,GXACluster.fetchCluster(scNVManager, accession, this, monitor));

		// Sanity check
	}

	public void fetchClusters () {
		new Thread(new FetchClusterThread()).start();
	}

	public void fetchDesign (final TaskMonitor monitor) {
		categories.set(1, GXADesign.fetchDesign(scNVManager, accession, this, monitor));

		// Sanity check
	}

	public void fetchDesign () {
		new Thread(new FetchDesignThread()).start();
	}

	// FIXME
	public void fetchIDF (final TaskMonitor monitor) {
	}

	public ZipInputStream getZipStream(String uri, TaskMonitor monitor) throws Exception {
		CloseableHttpClient httpclient = HttpClients.createDefault();
		HttpGet httpGet = new HttpGet(uri);
		CloseableHttpResponse response1 = httpclient.execute(httpGet);
		if (response1.getStatusLine().getStatusCode() != 200) {
			return null;
		}
		ZipInputStream stream = null;
		try {
			HttpEntity entity1 = response1.getEntity();
			stream = new ZipInputStream(entity1.getContent());
		} finally {
			response1.close();
		}
		return stream;
	}

	public String toString() {
		return gxaMetadata.toString();
	}

	class FetchClusterThread implements Runnable {
		@Override
		public void run() {
			categories.set(0, GXACluster.fetchCluster(scNVManager, accession, gxaExperiment, null));
		}
	}

	class FetchDesignThread implements Runnable {
		@Override
		public void run() {
			categories.set(1, GXADesign.fetchDesign(scNVManager, accession, gxaExperiment, null));
		}
	}
}
