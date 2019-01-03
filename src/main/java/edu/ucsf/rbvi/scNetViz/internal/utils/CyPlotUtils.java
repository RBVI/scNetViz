package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class CyPlotUtils {
	public static void createViolinPlot(ScNVManager manager, String names, String data, String title, 
	                                    String xlabel, String ylabel, String accession) {
		Map<String, Object> argMap = new HashMap<>();
		argMap.put("data", data);
		argMap.put("editor", "false");
		argMap.put("xlabel",xlabel);
		argMap.put("ylabel",ylabel);
		argMap.put("title",title);
		argMap.put("names",names);
		argMap.put("selectionString","scnetviz select accession=\""+accession+"\" genes=%s");
		manager.executeCommand("cyplot", "violin", argMap);
	}

	public static void createHeatMap(ScNVManager manager, String rowLabels, String columnLabels,
	                                 String data, String title, 
	                                 String xlabel, String ylabel, String accession) {
		Map<String, Object> argMap = new HashMap<>();
		argMap.put("rowLabels", rowLabels);
		argMap.put("columnLabels", columnLabels);
		argMap.put("data", data);
		argMap.put("editor", "false");
		argMap.put("xLabel",xlabel);
		argMap.put("yLabel",ylabel);
		argMap.put("title",title);
		argMap.put("selectionString","scnetviz select accession=\""+accession+"\" genes=%s");
		manager.executeCommand("cyplot", "heat", argMap);
	}

	// For Violin plots, we want to exclude NaNs, which means that our name array is different
	// for each trace.  So, we create a map with two strings, one for names and one for data
	public static String[] mapToDataAndNames(Map<String, double[]> dataMap, List<String> names) {
		StringBuilder dataBuilder = new StringBuilder();
		StringBuilder namesBuilder = new StringBuilder();
		dataBuilder.append("{");
		namesBuilder.append("{");
		for (String key: dataMap.keySet()) {

			dataBuilder.append("\""+key+"\":");
			dataBuilder.append("[");
			namesBuilder.append("\""+key+"\":");
			namesBuilder.append("[");
			double[] values = dataMap.get(key);
			for (int i = 0; i < values.length-1; i++) {
				if (!Double.isNaN(values[i])) {
					dataBuilder.append(String.valueOf(values[i])+",");
					namesBuilder.append(names.get(i)+",");
				}
			}
			dataBuilder.append(String.valueOf(values[values.length-1])+"],");
			namesBuilder.append(names.get(values.length-1)+"],");
		}
		String[] returnString = new String[2];
		returnString[0] = namesBuilder.substring(0, namesBuilder.length()-1)+"}";
		returnString[1] = dataBuilder.substring(0, dataBuilder.length()-1)+"}";
		return returnString;
	}

	public static String mapToData(Map<String, double[]> dataMap) {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		for (String key: dataMap.keySet()) {
			builder.append("\""+key+"\":");
			builder.append("[");
			double[] values = dataMap.get(key);
			for (int i = 0; i < values.length-1; i++) {
				if (Double.isNaN(values[i]))
					builder.append("\"NaN\",");
				else
					builder.append(String.valueOf(values[i])+",");
			}
			if (Double.isNaN(values[values.length-1]))
				builder.append("\"NaN\"],");
			else
				builder.append(String.valueOf(values[values.length-1])+"],");
		}
		return builder.substring(0, builder.length()-1)+"}";
	}

	public static String listToCSV(List<String> names) {
		StringBuilder builder = new StringBuilder();
		for (String s: names) {
			builder.append(s+",");
		}
		return builder.substring(0, builder.length()-1);
	}
}
