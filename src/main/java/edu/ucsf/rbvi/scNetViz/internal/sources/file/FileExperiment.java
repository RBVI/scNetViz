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

public class FileExperiment implements Experiment {
	final Logger logger;

	String accession = null;
	List<String[]> rowTable = null;
	List<String[]> colTable = null;
	MatrixMarket mtx = null;
	final List<Category> categories;
	// GXACluster fileCluster = null;
	// GXAIDF fileIDF = null;
	// GXADesign fileDesign = null;

	final ScNVManager scNVManager;
	final FileExperiment fileExperiment;
	final FileSource source;
	final FileMetadata fileMetadata;

	public FileExperiment (ScNVManager manager, FileSource source, FileMetadata metadata) {
		this.scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);
		this.fileExperiment = this;
		this.source = source;
		this.fileMetadata = metadata;
		categories = new ArrayList<Category>();
	}

	public Matrix getMatrix() { return mtx; }
	public String getAccession() { return accession; }

	public List<String[]> getColumnLabels() { return colTable; }
	public List<String[]> getRowLabels() { return rowTable; }

	// public GXACluster getClusters() { return fileCluster; }
	// public GXADesign getDesign() { return fileDesign; }

	public List<Category> getCategories() { return categories; }

	public void addCategory(Category c) { categories.add(c); }

	public Metadata getMetadata() { return fileMetadata; }

	public Source getSource() { return source; }

	public TableModel getTableModel() { return new FileExperimentTableModel(scNVManager, this); }

	public void readMTX (final TaskMonitor monitor) {
		// Get the URI
		File mtxFile = (File)fileMetadata.get(FileMetadata.FILE);

		try {
			FileInputStream inputStream = new FileInputStream(mtxFile);

			try {
				ZipInputStream zipStream = new ZipInputStream(inputStream);

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
				inputStream.close();
			}
		} catch (Exception e) {}
		scNVManager.addExperiment(accession, this);
	}

	public String toString() {
		return fileMetadata.toString();
	}

}
