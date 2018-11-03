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

import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.work.TaskMonitor;

public class CSVReader {
	public static String delimiter = null;

	public static List<String[]> readCSV(TaskMonitor taskMonitor, String name) throws IOException, FileNotFoundException {
		return readCSV(taskMonitor, new File(name));
	}

	public static List<String[]> readCSV(TaskMonitor taskMonitor, File name) throws IOException, FileNotFoundException {
		LogUtils.log(taskMonitor, TaskMonitor.Level.INFO, "Reading CSV file '"+name.toString()+"'");

		// Open the file
		BufferedReader input = new BufferedReader(new FileReader(name));
		return readCSV(taskMonitor, input);
	}

	public static List<String[]> readCSV(TaskMonitor taskMonitor, InputStream stream, String name) throws IOException {
		LogUtils.log(taskMonitor, TaskMonitor.Level.INFO, "Reading CSV file '"+name.toString()+"'");
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

		LogUtils.log(taskMonitor, TaskMonitor.Level.INFO, "Found "+rowList.size()+" rows with "+rowList.get(0).length+" labels each");
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
				LogUtils.log(taskMonitor, TaskMonitor.Level.ERROR, 
				    "Error return from '"+fetchString+"': "+response1.getStatusLine());
				return null;
			}
			HttpEntity entity1 = response1.getEntity();

			try {
				InputStream inputStream = entity1.getContent();
				input = readCSV(taskMonitor, inputStream, accession);
				inputStream.close();
			} catch (Exception e) {
				LogUtils.log(taskMonitor, TaskMonitor.Level.ERROR, "Error reading from '"+fetchString+"': "+e.getMessage());
			} finally {
				response1.close();
			}
		} catch (Exception e) {
			LogUtils.log(taskMonitor, TaskMonitor.Level.ERROR, 
			    "Error attempting to fetch '"+fetchString+"': "+e.getMessage());
		}
		return input;
	}

	private static String[] readRow(BufferedReader input) throws IOException {
		String row = input.readLine();
		// System.out.println("Row: "+row);
		if (row == null) return null;
		return smartSplit(row);
	}

	public static String[] smartSplit(String input) {
		String[] splitString;
		if (delimiter != null)
			splitString = input.split(delimiter, -1);
		else {
			delimiter = "\t";
			splitString = input.split(delimiter, -1);
			if (splitString.length == 1) {
				delimiter = ",";
				splitString = input.split(delimiter, -1);
				if (splitString.length == 1) {
					delimiter = null;
					throw new RuntimeException("Only tabs and commas are supported delimiters");
				}
			}
		}
		return splitString;
	}

}
