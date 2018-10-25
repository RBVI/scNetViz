package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;
import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.model.MTXManager;

public class CSVReader {
	public static String delimiter = null;
	static Logger logger = Logger.getLogger(CyUserLog.NAME);

	public static List<String[]> readCSV(TaskMonitor taskMonitor, String name) throws IOException, FileNotFoundException {
		return readCSV(taskMonitor, new File(name));
	}

	public static List<String[]> readCSV(TaskMonitor taskMonitor, File name) throws IOException, FileNotFoundException {
		log(taskMonitor, TaskMonitor.Level.INFO, "Reading CSV file '"+name.toString()+"'");

		// Open the file
		BufferedReader input = new BufferedReader(new FileReader(name));
		return readCSV(taskMonitor, input);
	}

	public static List<String[]> readCSV(TaskMonitor taskMonitor, InputStream stream, String name) throws IOException {
		log(taskMonitor, TaskMonitor.Level.INFO, "Reading CSV file '"+name.toString()+"'");
		BufferedReader input = new BufferedReader(new InputStreamReader(stream));
		return readCSV(taskMonitor, input);
	}

	public static List<String[]> readCSV(TaskMonitor taskMonitor, BufferedReader input) throws IOException, FileNotFoundException {

		// Read each row
		List<String[]> rowList = new ArrayList<>();
		do {
			String[] row = readRow(input);
			if (row == null)
				break;
			rowList.add(row);
		} while (true);

		log(taskMonitor, TaskMonitor.Level.INFO, "Found "+rowList.size()+" rows with "+rowList.get(0).length+" labels each");
		return rowList;
	}

	public static List<String[]> readCSVFromHTTP(TaskMonitor taskMonitor, 
	                                             String URI, String accession) {
		List<String[]> input = null;
		String fetchString = String.format(URI, accession);
		try {
			CloseableHttpClient httpclient = HttpClients.createDefault();
			HttpGet httpGet = new HttpGet(fetchString);
			CloseableHttpResponse response1 = httpclient.execute(httpGet);
			if (response1.getStatusLine().getStatusCode() != 200) {
				log(taskMonitor, TaskMonitor.Level.ERROR, 
				    "Error return from '"+fetchString+"': "+response1.getStatusLine());
				return null;
			}
			HttpEntity entity1 = response1.getEntity();

			try {
				InputStream inputStream = entity1.getContent();
				input = readCSV(taskMonitor, inputStream, accession);
				inputStream.close();
			} catch (Exception e) {
				log(taskMonitor, TaskMonitor.Level.ERROR, "Error reading from '"+fetchString+"': "+e.getMessage());
			} finally {
				response1.close();
			}
		} catch (Exception e) {
			log(taskMonitor, TaskMonitor.Level.ERROR, 
			    "Error attempting to fetch '"+fetchString+"': "+e.getMessage());
		}
		return input;
	}

	private static String[] readRow(BufferedReader input) throws IOException {
		String row = input.readLine();
		// System.out.println("Row: "+row);
		if (row == null) return null;
		String[] columns;
		if (delimiter != null)
			columns = row.split(delimiter, -1);
		else {
			delimiter = "\t";
			columns = row.split(delimiter, -1);
			if (columns.length == 1) {
				delimiter = ",";
				columns = row.split(delimiter, -1);
				if (columns.length == 1) {
					delimiter = null;
					throw new RuntimeException("Only tabs and commas are supported column delimiters");
				}
			}
		}
		return columns;
	}

	private static void log(TaskMonitor taskMonitor, TaskMonitor.Level level, String message) {
		if (taskMonitor != null) {
			taskMonitor.showMessage(level, message);
			return;
		}
		switch (level) {
			case ERROR:
				logger.error(message);
				break;
			case INFO:
				logger.info(message);
				break;
			case WARN:
				logger.warn(message);
				break;
		}
	}
}
