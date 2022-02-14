package edu.ucsf.rbvi.scNetViz.internal.sources.file;

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
import java.util.zip.GZIPInputStream;

import org.apache.log4j.Logger;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;

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
import edu.ucsf.rbvi.scNetViz.internal.utils.FileUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

public class FileExperiment implements Experiment {
	final Logger logger;

	String accession = null;
	List<String[]> rowTable = null;
	List<String[]> colTable = null;
	MatrixMarket mtx = null;
	final List<Category> categories;
	double[][] tSNE;
	String plotType = null;
	// GXACluster fileCluster = null;
	// GXAIDF fileIDF = null;
	// GXADesign fileDesign = null;

	final ScNVManager scNVManager;
	final FileExperiment fileExperiment;
	final FileSource source;
	final FileMetadata fileMetadata;
	DifferentialExpression diffExp = null;
	FileExperimentTableModel tableModel = null;

	int rowIndexKey = 1;
	int columnIndexKey = 1;

	public static String SERVICES_URI = "https://webservices.rbvi.ucsf.edu/scnetviz/api/v2/save/File/%s";

	public FileExperiment (ScNVManager manager, FileSource source, FileMetadata metadata) {
		this.scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);
		this.fileExperiment = this;
		this.source = source;
		this.fileMetadata = metadata;
		this.accession = metadata.get(Metadata.ACCESSION).toString();
		categories = new ArrayList<Category>();
	}

	public Matrix getMatrix() { return mtx; }
	public String getAccession() { return accession; }

	public List<String[]> getColumnLabels() { return colTable; }
	public List<String[]> getRowLabels() { return rowTable; }

	// public GXACluster getClusters() { return fileCluster; }
	// public GXADesign getDesign() { return fileDesign; }

	public List<Category> getCategories() { return categories; }

	public Category getCategory(String categoryName) { 
		for (Category cat: categories) {
			if (cat.toString().equals(categoryName))
				return cat;
		}
		return null;
	}

	public void addCategory(Category c) { categories.add(c); }

	public Category getDefaultCategory() { 
		if (categories.size() > 0)
			return categories.get(0); 
		return null;
	}

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

	public Metadata getMetadata() { return fileMetadata; }

	public Source getSource() { return source; }

	public String getSpecies() { return (String)fileMetadata.get(Metadata.SPECIES); }

	public FileExperimentTableModel getTableModel() { 
		if (tableModel == null)
			tableModel = new FileExperimentTableModel(scNVManager, this); 
		return tableModel;
	}

	public DifferentialExpression getDiffExp() { return diffExp; }
	public void setDiffExp(DifferentialExpression de) { diffExp = de; }

	public void readMTX (final TaskMonitor monitor, boolean skipFirst) {
		// Initialize
		mtx = null;
		colTable = null;
		rowTable = null;

		// Get the URI
		File mtxFile = (File)fileMetadata.get(FileMetadata.FILE);

		// Three possibilities:
		// 1) Directory
		// 2) ZIP file
		// 3) Tar.gz file

		try {
			if (mtxFile.isDirectory()) {
				// System.out.println("Directory");
				for (File f: mtxFile.listFiles()) {
					readFile(monitor, f, skipFirst);
				}
			} else if (FileUtils.isZip(mtxFile.getName())) {
				// System.out.println("Zip file");
				ZipInputStream zipStream = FileUtils.getZipInputStream(mtxFile);
				ZipEntry entry;
				while ((entry = zipStream.getNextEntry()) != null) {
					readFile(monitor, entry.getName(), zipStream, skipFirst);
					zipStream.closeEntry();
				}
				zipStream.close();
			} else if (FileUtils.isTar(mtxFile.getName())) {
				// System.out.println("Tar file");
				TarArchiveInputStream tarStream = FileUtils.getTarInputStream(mtxFile);
				TarArchiveEntry entry;
				while ((entry = tarStream.getNextTarEntry()) != null) {
					// System.out.println("Tar entry: "+entry.getName());
					if (FileUtils.isGzip(entry.getName())) {
						InputStream stream = FileUtils.getGzipStream(tarStream);
						readFile(monitor, entry.getName(), stream, skipFirst);
					} else {
						readFile(monitor, entry.getName(), tarStream, skipFirst);
					}
				}
				tarStream.close();
			} else {
				readFile(monitor, mtxFile, skipFirst);
			}
    } catch (FileNotFoundException e) {
      monitor.showMessage(TaskMonitor.Level.ERROR, "No such file: '"+mtxFile.getName()+"'");
      return;
		} catch(IOException e) {
      monitor.showMessage(TaskMonitor.Level.ERROR, "Error reading file: '"+mtxFile.getName()+"': "+e.getMessage());
      e.printStackTrace();
      return;
    }

		scNVManager.addExperiment(accession, this);

    // Cache the file on the server
		new Thread(new CacheExperimentThread(accession)).start();

		System.out.println("mtx has "+mtx.getNRows()+" rows and "+mtx.getNCols()+" columns");
	}

	private void readFile(TaskMonitor monitor, File f, boolean skipFirst) throws IOException {
		InputStream stream = new FileInputStream(f);
		// if (FileUtils.isGzip(f.getName())) {
		// 	stream = FileUtils.getGzipStream(stream);
		// }
		readFile(monitor, f.getName(), stream, skipFirst);
	}

	private void readFile(TaskMonitor monitor, String name, InputStream stream, 
	                      boolean skipFirst) throws IOException {
		if (isColumnFile(name)) {
			// System.out.println("Reading columns from "+name);
			colTable = CSVReader.readCSV(monitor, stream, name);
			// System.out.println("colTable has "+colTable.size()+" columns");

			if (skipFirst) {
				if (colTable.size() > 1)
					colTable.remove(0);
				// See if the first line is a header
				// FileUtils.skipHeader(colTable);
			}

			columnIndexKey = getColumnIndex(colTable, name);
			if (mtx != null) {
				mtx.setColumnTable(colTable, columnIndexKey);
			}
		} else if (isRowFile(name)) {
			rowTable = CSVReader.readCSV(monitor, stream, name);
			// System.out.println("rowTable has "+rowTable.size()+" rows");
			if (skipFirst) {
				if (rowTable.size() > 1)
					rowTable.remove(0);
				// See if the first line is a header
				// FileUtils.skipHeader(rowTable);
			}

			rowIndexKey = getRowIndex(rowTable, name);
			if (mtx != null) {
				mtx.setRowTable(rowTable, rowIndexKey);
			}
		} if (isMtxFile(name)) {
			mtx = new MatrixMarket(scNVManager, null, null);
			if (rowTable != null)
				mtx.setRowTable(rowTable, rowIndexKey);
			if (colTable != null)
				mtx.setColumnTable(colTable, columnIndexKey);

			mtx.readMTX(monitor, stream, name);
		}
	}

	public String toString() {
		return getAccession();
	}

	public String toHTML() {
		return fileMetadata.toHTML();
	}

	public String toJSON() {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		builder.append("\"source\": \""+getSource().toString()+"\",\n");
		builder.append("\"source name\": \""+getSource().getName()+"\",\n");
		builder.append("\"metadata\": "+fileMetadata.toJSON()+",\n");
		builder.append("\"rows\": "+getMatrix().getNRows()+",\n");
		builder.append("\"columns\": "+getMatrix().getNCols()+",\n");
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

	private boolean isColumnFile(String fileName) {
		String name = FileUtils.baseName(fileName);
		if (name.endsWith(".mtx_cols") || name.contains("colLabels") || name.startsWith("barcodes")) {
			return true;
		}
		return false;
	}

	private boolean isRowFile(String fileName) {
		String name = FileUtils.baseName(fileName);
		if (name.endsWith(".mtx_rows") || name.contains("rowLabels") || name.startsWith("features")) {
			return true;
		}
		return false;
	}

	private boolean isMtxFile(String fileName) {
		if (fileName.endsWith(".mtx") || fileName.endsWith(".mtx.gz"))
			return true;
		return false;
	}

	private int getColumnIndex(List<String[]> table, String name) {
		// System.out.println("getColumnIndex of "+name);
		// System.out.println("getColumnIndex table size = "+table.size());
		// System.out.println("getColumnIndex table length = "+table.get(0).length);
		if (table.size() > 0 && table.get(0).length == 1 && name.contains("barcodes"))
			return 0;
		return 1;
	}

	private int getRowIndex(List<String[]> table, String name) {
		// System.out.println("Row Table size = "+table.size());
		if (table.size() > 0 && name.contains("features"))
			return 0;
		return 1;
	}

	class CacheExperimentThread implements Runnable {
    final String accession;
    public CacheExperimentThread(String acc) { this.accession = acc; }
		@Override
		public void run() {
      String postString = String.format(SERVICES_URI, accession);
      try {
        // Create the cache file
        System.out.println("Creating cache file");
        mtx.createCache("File", accession);

        // Wait until the cache is available
        while (!mtx.hasCache()) {
          Thread.sleep(1000);
          System.out.println("mtx.hasCache = "+mtx.hasCache());
        }

        System.out.println("Sending file to server");
        // OK, now send the file to the server
        File expFile = mtx.getMatrixCache();
        HTTPUtils.postFile(postString, expFile, null);
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

}
