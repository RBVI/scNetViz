package edu.ucsf.rbvi.scNetViz.internal.sources.file;

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
import edu.ucsf.rbvi.scNetViz.internal.utils.FileUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

public class FileExperiment implements Experiment {
	final Logger logger;

	String accession = null;
	List<String[]> rowTable = null;
	List<String[]> colTable = null;
	MatrixMarket mtx = null;
	final List<Category> categories;
	double[][] tSNE;
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
				System.out.println("Directory");
				for (File f: mtxFile.listFiles()) {
					readFile(monitor, f, skipFirst);
				}
			} else if (FileUtils.isZip(mtxFile.getName())) {
				System.out.println("Zip file");
				ZipInputStream zipStream = FileUtils.getZipInputStream(mtxFile);
				ZipEntry entry;
				while ((entry = zipStream.getNextEntry()) != null) {
					readFile(monitor, entry.getName(), zipStream, skipFirst);
					zipStream.closeEntry();
				}
				zipStream.close();
			} else if (FileUtils.isTar(mtxFile.getName())) {
				System.out.println("Tar file");
				TarArchiveInputStream tarStream = FileUtils.getTarInputStream(mtxFile);
				TarArchiveEntry entry;
				while ((entry = tarStream.getNextTarEntry()) != null) {
					System.out.println("Tar entry: "+entry.getName());
					if (FileUtils.isGzip(entry.getName())) {
						InputStream stream = FileUtils.getGzipStream(tarStream);
						readFile(monitor, entry.getName(), stream, skipFirst);
					} else {
						readFile(monitor, entry.getName(), tarStream, skipFirst);
					}
				}
				tarStream.close();
			}
		} catch(IOException e) {}

		scNVManager.addExperiment(accession, this);
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
		builder.append("source: '"+getSource().toString()+"',");
		builder.append("accession: '"+getMetadata().get(Metadata.ACCESSION).toString()+"',");
		builder.append("species: '"+getSpecies().toString()+"',");
		builder.append("description: '"+getMetadata().get(Metadata.DESCRIPTION).toString()+"',");
		builder.append("rows: '"+getMatrix().getNRows()+"',");
		builder.append("columns: '"+getMatrix().getNCols()+"',");
		List<Category> categories = getCategories();
		builder.append("categories: [");
		for (Category cat: categories) {
			builder.append(cat.toJSON()+",");
		}
		return builder.substring(0, builder.length()-1)+"]}";
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
		System.out.println("getColumnIndex of "+name);
		System.out.println("getColumnIndex table size = "+table.size());
		System.out.println("getColumnIndex table length = "+table.get(0).length);
		if (table.size() > 0 && table.get(0).length == 1 && name.contains("barcodes"))
			return 0;
		return 1;
	}

	private int getRowIndex(List<String[]> table, String name) {
		System.out.println("Row Table size = "+table.size());
		if (table.size() > 0 && name.contains("features"))
			return 0;
		return 1;
	}

}
