package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
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

import org.apache.http.Header;
import org.apache.http.HttpEntity;
import org.apache.http.client.config.RequestConfig;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
import org.apache.log4j.Logger;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskMonitor;

import org.json.simple.JSONObject;

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

public class GXAExperiment implements Experiment {
	public static String RESULTS_URL = "https://www.ebi.ac.uk/gxa/sc/experiments/%s/Results";
	public static String GXA_MTX_URI = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download/zip?fileType=normalised";
	public static String SERVICES_URI = "http://webservices.rbvi.ucsf.edu/scnetviz/api/v2/fetch/GXA/%s";
	final Logger logger;

	String accession = null;
	List<String[]> rowTable = null;
	List<String[]> colTable = null;
	MatrixMarket mtx = null;
	GXAMetadata gxaMetadata = null;
	List<Category> categories;
	double[][] tSNE;
	String plotType = null;

	final ScNVManager scNVManager;
	final GXAExperiment gxaExperiment;
	final GXASource source;
	GXAExperimentTableModel tableModel = null;

	DifferentialExpression diffExp = null;

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

	public Category getCategory(String categoryName) { 
		for (Category cat: categories) {
			if (cat.toString().equals(categoryName))
				return cat;
		}
		return null;
	}

	public void addCategory(Category c) { categories.add(c); }

	public Metadata getMetadata() { return gxaMetadata; }

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
			tableModel = new GXAExperimentTableModel(scNVManager, this);
		return tableModel;
	}

	public void fetchMTX (final TaskMonitor monitor) {
		// Get the URI
		String fetchString = String.format(GXA_MTX_URI, accession);

		// IMPORTANT!  We need to do a gc here to avoid any pauses later
		Runtime.getRuntime().gc();

    // Spawn off a thread to cache this accession on the server
    cacheAccession();

		// See how many cells we have.  For really large matrices, we want to "peak" first and
		// preallocate our arrays.
		int assays = ((Long)gxaMetadata.get(GXAMetadata.ASSAYS)).intValue();

		if (assays > 40000) {
			fetchAccession(monitor, fetchString, true);
		}

		fetchAccession(monitor, fetchString, false);

		// if (mtx != null)
		// 	mtx.createCache(source.getName(), accession);

		scNVManager.addExperiment(accession, this);

	}

	private void fetchAccession(TaskMonitor monitor, String fetchString, boolean peak) {
		// System.out.println("fetchAccession");
		try {
			CloseableHttpClient httpclient = HttpClients.createDefault();
			RequestConfig config = RequestConfig.custom()
							.setSocketTimeout(60*1000)
							.build();
			HttpGet httpGet = new HttpGet(fetchString);
			httpGet.setConfig(config);
			CloseableHttpResponse response1 = httpclient.execute(httpGet);
			if (response1.getStatusLine().getStatusCode() != 200) {
				return;
			}
			/*
			Header[] headers = response1.getAllHeaders();
			for (Header header: headers) {
				System.out.println("Key: "+header.getName()+" ,Value: "+header.getValue());
			}
			*/
			HttpEntity entity1 = response1.getEntity();

			try {
				ZipInputStream zipStream = new ZipInputStream(entity1.getContent());

				ZipEntry entry;
				while ((entry = zipStream.getNextEntry()) != null) {
					String name = entry.getName();

					if (peak) {
						// Use peak to pre-allocate all of the necessary memory
						if (name.endsWith(".mtx")) {
							mtx = new MatrixMarket(scNVManager, null, null);
							mtx.peak(monitor, zipStream, name);
							zipStream.closeEntry();
							zipStream.close();
							return;
						} else {
							continue;
						}
					}

					if (name.endsWith(".mtx_cols")) {
						colTable = CSVReader.readCSV(monitor, zipStream, name);
						if (mtx != null) 
							mtx.setColumnTable(colTable,0);
					} else if (name.endsWith(".mtx_rows")) {
						rowTable = CSVReader.readCSV(monitor, zipStream, name);
						if (mtx != null) 
							mtx.setRowTable(rowTable, 0);
					} else if (name.endsWith(".mtx")) {
						if (mtx == null)
							mtx = new MatrixMarket(scNVManager, null, null);
						mtx.setRowTable(rowTable);
						mtx.setColumnTable(colTable);
						mtx.readMTX(monitor, zipStream, name);
						/*
						{
							FileOutputStream outStream = new FileOutputStream("/tmp/"+name);
							byte[] buffer = new byte[8192];
							int count = 0;
							while ((count = zipStream.read(buffer)) >= 0) {
								outStream.write(buffer, 0, count);
							}
							outStream.close();
						}
						*/
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
	}


	public void fetchClusters (final TaskMonitor monitor) {
		categories.set(0,GXACluster.fetchCluster(scNVManager, accession, this, monitor));

		// Sanity check
	}

  private void cacheAccession() {
		new Thread(new CacheExperimentThread()).start();
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

	public String toHTML() {
		return gxaMetadata.toHTML();
	}

	public String toString() {
		return getAccession();
	}

	public String toJSON() {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		builder.append("\"source name\": \""+getSource().getName()+"\",\n");
		builder.append("\"source\": \""+getSource().toString()+"\",\n");
		builder.append("\"metadata\": "+gxaMetadata.toJSON()+",\n");
		builder.append("\"rows\": "+getMatrix().getNRows()+",\n");
		builder.append("\"columns\": "+getMatrix().getNCols()+",\n");
		List<Category> categories = getCategories();
		builder.append("\"categories\": [");
		int nCat = categories.size();
		for (Category cat: categories) {
			builder.append(cat.toJSON());
			if (nCat-- > 1) builder.append(",\n");
		}
		builder.append("]");
		if (diffExp != null) {
			builder.append(",\""+ScNVManager.DIFFEXP+"\":");
			builder.append(diffExp.toJSON()+"\n");
		}
		builder.append("}");
		return builder.toString();
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

	public void loadFromSession(Map<String, File> fileMap) throws IOException {
		// Find our files
		String expPrefix = URLEncoder.encode(source.getName()+"."+accession);
		if (!fileMap.containsKey(expPrefix+".mtx"))
			throw new FileNotFoundException("File '"+expPrefix+".mtx' doesn't exist");

		File mtxFile = fileMap.get(expPrefix+".mtx");

		if (!fileMap.containsKey(expPrefix+".mtx_rows"))
			throw new FileNotFoundException("File '"+expPrefix+".mtx_rows' doesn't exist");

		File rowFile = fileMap.get(expPrefix+".mtx_rows");

		if (!fileMap.containsKey(expPrefix+".mtx_cols"))
			throw new FileNotFoundException("File '"+expPrefix+".mtx_cols' doesn't exist");

		File colFile = fileMap.get(expPrefix+".mtx_cols");

		// Read in the MatrixMarket file
		System.out.println("Reading MTX file");
		mtx = new MatrixMarket(scNVManager);
		mtx.readMTX(null, fileMap.get(expPrefix+".mtx"));
	
		System.out.println("Reading column header file");
		// Read in our row and column headers
		colTable = CSVReader.readCSV(null, colFile);
		rowTable = CSVReader.readCSV(null, rowFile);
		mtx.setRowTable(rowTable);
		mtx.setColumnTable(colTable, 0);
		scNVManager.addExperiment(accession, this);
	}

	public Category loadCategoryFromSession(JSONObject jsonCategory, Map<String, File> fileMap) throws IOException {
		// See which category it is
		String name = (String)jsonCategory.get("name");
		System.out.println("Loading category: '"+name+"' from session");
		String catSource = (String)jsonCategory.get("source name");
		System.out.println("Source = "+catSource);
		String catPrefix = URLEncoder.encode(source.getName()+"."+accession+"."+catSource+"."+name);
		String fileName = catPrefix+".csv";
		System.out.println("fileName = "+fileName);
		if (!fileMap.containsKey(fileName))
			throw new FileNotFoundException("File '"+fileName+"' doesn't exist");

		if (name.equals("Cluster")) {
			System.out.println("Reading cluster");
			GXACluster cluster = GXACluster.readCluster(scNVManager, this, fileMap.get(fileName), jsonCategory);
			categories.set(0, cluster);
			return cluster;
		} else if (name.equals("Design/Factors")) {
			System.out.println("Reading design/factors");
			GXADesign design = GXADesign.readDesign(scNVManager, this, fileMap.get(fileName), jsonCategory);
			categories.set(1, design);
			return design;
		}
		return null;
	}

	public DifferentialExpression loadDiffExpFromSession(JSONObject jsonDiffExp, Map<String, File> fileMap) throws IOException {
		// Get the file
		String deFileName = URLEncoder.encode(source.getName()+"."+accession+".diffExp")+".csv";
		if (fileMap.containsKey(deFileName)) {
			File deFile = fileMap.get(deFileName);
			try {
				diffExp = new DifferentialExpression(scNVManager, this, jsonDiffExp, deFile);
			} catch (Exception e) {
				logger.error("Unable to read differential expression data for "+accession+" in session: "+e.toString());
				e.printStackTrace();
				return null;
			}
			return diffExp;
		}
		return null;
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

  // Force the server to cache the experiment for us
	class CacheExperimentThread implements Runnable {
		@Override
		public void run() {
      String fetchString = String.format(SERVICES_URI, accession);
      try {
        CloseableHttpClient httpclient = HttpClients.createDefault();
        RequestConfig config = RequestConfig.custom()
                .setSocketTimeout(60*1000)
                .build();
        HttpGet httpGet = new HttpGet(fetchString);
        httpGet.setConfig(config);
        CloseableHttpResponse response1 = httpclient.execute(httpGet);
        if (response1.getStatusLine().getStatusCode() == 200) {
          return;
        }
			} catch (Exception e) {
				e.printStackTrace();
      }
		}
	}
}
