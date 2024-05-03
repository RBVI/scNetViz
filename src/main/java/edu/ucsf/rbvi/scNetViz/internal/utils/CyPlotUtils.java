package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class CyPlotUtils {
	public static void createViolinPlot(ScNVManager manager, String names, String data, 
	                                    String groups, String title, 
	                                    String xlabel, String ylabel, String accession, 
	                                    boolean showAll) {
		createViolinPlot(manager, names, data, groups, title, xlabel, ylabel, accession, showAll, 0.0);
	}

	public static void createViolinPlot(ScNVManager manager, String names, String data, 
	                                    String groups, String title, 
	                                    String xlabel, String ylabel, String accession, 
	                                    boolean showAll, double jitter) {
		// System.out.println("createViolinPlot");
		Map<String, Object> argMap = new HashMap<>();
		argMap.put("data", data);
		argMap.put("editor", "false");
		argMap.put("xLabel",xlabel);
		argMap.put("yLabel",ylabel);
		argMap.put("title",title);
		argMap.put("names",names);
		argMap.put("groups",groups);
		argMap.put("id",accession);
		if (showAll)
			argMap.put("showAll","true");
		else
			argMap.put("showAll","false");
		argMap.put("jitter",jitter);
		argMap.put("selectionString","scnetviz select accession=\""+accession+"\" genes=%s");
		manager.executeCommand("cyplot", "violin", argMap);
	}

	public static void createScatterPlot(ScNVManager manager, String names, String xValues, String yValues,
	                                     String zValues,
	                                     String title, 
	                                     String xlabel, String ylabel, String accession) {
		Map<String, Object> argMap = new HashMap<>();
		argMap.put("xValues", xValues);
		argMap.put("yValues", yValues);
		argMap.put("editor", "false");
		argMap.put("xLabel",xlabel);
		argMap.put("yLabel",ylabel);
		if (zValues != null) {
			argMap.put("zValues",zValues);
			argMap.put("colorscale","Viridis");
			argMap.put("editor","false");
			argMap.put("scaleLabel", "log(TPM)");
		} else {
			argMap.put("editor","false");
			argMap.put("scaleLabel", "TPM");
		}
		argMap.put("title",title);
		argMap.put("names",names);
		argMap.put("id",accession);
		argMap.put("selectionString","scnetviz select accession=\""+accession+"\" cells=%s");
		manager.executeCommand("cyplot", "scatter", argMap);
	}

	public static void createHeatMap(ScNVManager manager, String rowLabels, String columnLabels,
	                                 String data, String title, 
	                                 String xlabel, String ylabel, String accession, boolean posOnly) {
		Map<String, Object> argMap = new HashMap<>();
		argMap.put("rowLabels", rowLabels);
		argMap.put("columnLabels", columnLabels);
		argMap.put("data", data);
		argMap.put("editor", "false");
		argMap.put("xLabel",xlabel);
		argMap.put("yLabel",ylabel);
		argMap.put("title",title);
		argMap.put("id",accession);
		if (posOnly)
			argMap.put("palette", "Viridis Sequential Viridis perceptually balanced palette");
		argMap.put("selectionString","scnetviz select accession=\""+accession+"\" genes=%s");
		manager.executeCommand("cyplot", "heat", argMap);
	}

	// For Violin plots, we want to exclude NaNs, which means that our name array is different
	// for each trace.  So, we create a map with two strings, one for names and one for data
	public static String[] mapToDataAndNames(Map<String, double[]> dataMap, List<String> names, 
	                                         List<String> columns) {
		StringBuilder dataBuilder = new StringBuilder();
		StringBuilder namesBuilder = new StringBuilder();
		dataBuilder.append("{");
		namesBuilder.append("{");
		for (String key: columns) {
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
			if (!Double.isNaN(values[values.length-1])) {
				dataBuilder.append(String.valueOf(values[values.length-1])+"],");
				namesBuilder.append(names.get(values.length-1)+"],");
			} else {
				dataBuilder.append("],");
				namesBuilder.append("],");
			}
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

	public static String listToJSON(List<String[]> names, int hdr) {
		StringBuilder builder = new StringBuilder();
		builder.append("[");
		for (String[] s: names) {
			builder.append("\""+s[hdr]+"\",");
		}
		return builder.substring(0, builder.length()-1)+"]";
	}

	public static String coordinatesToJSON(Map<String,double[]> coords, int index) {
		StringBuilder builder = new StringBuilder();
		builder.append("[");

		// This assumes that the coords map is an ordered map (i.e. a LinkedHashMap)
		for (String cell: coords.keySet()) {
			double[] coord = coords.get(cell);
			builder.append(coord[index]+",");
		}
		builder.setCharAt(builder.length()-1, ']');
		return builder.toString();
	}

	/*
	public static String coordinatesToJSON(double[][] coords, int index) {
		StringBuilder builder = new StringBuilder();
		builder.append("[");
		for (int cell = 0; cell < coords.length-1; cell++) {
			builder.append(coords[cell][index]+",");
		}
		builder.append(coords[coords.length-1][index]+"]");
		return builder.toString();
	}
	*/

	public static String valuesToJSON(DoubleMatrix matrix, int row, boolean log) {
		StringBuilder builder = new StringBuilder();
		builder.append("[");
		for (int column = 0; column < matrix.getNCols()-1; column++) {
			builder.append(getNonNaNValue(matrix.getDoubleValue(row, column),0.0, log)+",");
		}
		builder.append(getNonNaNValue(matrix.getDoubleValue(row, matrix.getNCols()-1),0.0, log)+"]");
		return builder.toString();
	}

	public static String listToMap(Map<Object, List<String>> map, List<String[]> list, int hdr) {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		for (Object trace: map.keySet()) {
			builder.append("\""+trace.toString()+"\":[");
			for (String cell: map.get(trace)) {
				builder.append("\""+cell+"\",");
			}
			builder.setCharAt(builder.length()-1, ']');
			builder.append(",");
		}
		builder.setCharAt(builder.length()-1, '}');
		return builder.toString();
	}

	public static String coordsToMap(Map<Object, List<String>> map, Map<String,double[]> coords, int index) {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		for (Object trace: map.keySet()) {
			builder.append("\""+trace.toString()+"\":[");
			for (String cell: map.get(trace)) {
				if (coords.containsKey(cell)) {
					builder.append(coords.get(cell)[index]+",");
				}
			}
			builder.setCharAt(builder.length()-1, ']');
			builder.append(",");
		}
		builder.setCharAt(builder.length()-1, '}');
		return builder.toString();
	}

	/*
	public static String coordsToMap(Map<Object, List<Integer>> map, double[][] coords, int index) {
		StringBuilder builder = new StringBuilder();
		builder.append("{");
		for (Object trace: map.keySet()) {
			builder.append("\""+trace.toString()+"\":[");
			for (Integer i: map.get(trace)) {
				builder.append(coords[i][index]+",");
			}
			builder.setCharAt(builder.length()-1, ']');
			builder.append(",");
		}
		builder.setCharAt(builder.length()-1, '}');
		return builder.toString();
	}
	*/

	public static double getNonNaNValue(double v, double missing, boolean log) {
		if (!Double.isNaN(v))
			if (log)
				return Math.log(v);
			else
				return v;
		return missing;
	}
}
