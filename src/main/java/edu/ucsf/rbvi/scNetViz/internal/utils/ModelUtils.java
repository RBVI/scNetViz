package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.awt.Color;
import java.awt.Paint;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import org.cytoscape.model.CyIdentifiable;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.CyTableUtil;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.util.color.Palette;
import org.cytoscape.util.color.PaletteProvider;
import org.cytoscape.util.color.PaletteProviderManager;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualMappingManager;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.view.vizmap.VisualStyleFactory;
import org.cytoscape.view.vizmap.mappings.BoundaryRangeValues;
import org.cytoscape.view.vizmap.mappings.ContinuousMapping;
import org.cytoscape.work.json.JSONResult;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class ModelUtils {
	public static final String EXPERIMENT_ACCESSION = "Experiment Accession";
	public static final String EXPERIMENT_SOURCE = "Experiment Source";
	public static final String CATEGORY = "Category";
	public static final String CATEGORY_NAMES = "Category Names";
	public static final String CATEGORY_ROW = "Category Row";
	public static final String CATEGORY_SOURCE = "Category Source";

	public static final String POSITIVE_ONLY = "positiveOnly";
	public static final String NEGATIVE_ONLY = "negativeOnly";
	public static final String SELECTED_ONLY = "selectedOnly";
	public static final String ENTIRE_NETWORK = "entireNetwork";

	public static final String NAMESPACE = "scNetViz";
;
	public static CyNetwork getNetworkFromJSON(ScNVManager manager, JSONResult res) {
		try {
			JSONObject json = (JSONObject) new JSONParser().parse(res.getJSON());
			long networkSUID = 0;
			if (json.containsKey("SUID")) {
				networkSUID = ((Number)json.get("SUID")).longValue();
			} else if (json.containsKey("network")) {
				networkSUID = ((Number)json.get("network")).longValue();
			}

			CyNetwork network = manager.getService(CyNetworkManager.class).getNetwork(networkSUID);
			return network;
		} catch (Exception e) { return null; }
	}

	public static void rename(CyNetwork network, CyIdentifiable id, String name) {
		network.getRow(id).set(CyNetwork.NAME, name);
	}

	public static String getName(CyNetwork network, CyIdentifiable id) {
		return network.getRow(id).get(CyNetwork.NAME, String.class);
	}

	public static void createDEColumns(ScNVManager manager, CyNetwork network, DifferentialExpression diffExp, 
	                                   String label) {
		CyTable nodeTable = network.getDefaultNodeTable();
		createColumnIfNeeded(nodeTable, NAMESPACE, label+" MTC", Double.class);
		createColumnIfNeeded(nodeTable, NAMESPACE, label+" Min.pct", Double.class);
		createColumnIfNeeded(nodeTable, NAMESPACE, label+" MDTC", Double.class);
		createColumnIfNeeded(nodeTable, NAMESPACE, label+" log2FC", Double.class);
		createColumnIfNeeded(nodeTable, NAMESPACE, label+" pValue", Double.class);
		createColumnIfNeeded(nodeTable, NAMESPACE, label+" FDR", Double.class);
		createColumnIfNeeded(nodeTable, NAMESPACE, label+" Rank", Integer.class);
	}

	public static void updateDEData(ScNVManager manager, CyNetwork network, 
	                                List<String> geneList, DifferentialExpression diffExp, 
	                                String label, List<String> rankedGenes) {
		Map<String, CyRow> queryMap = createQueryMap(network);
		for (String gene: geneList) {
			CyRow row = queryMap.get(gene);
			if (row == null) continue;

			updateValues(network, row, diffExp, gene, label+" MTC");
			updateValues(network, row, diffExp, gene, label+" Min.pct");
			updateValues(network, row, diffExp, gene, label+" MDTC");
			updateValues(network, row, diffExp, gene, label+" log2FC");
			updateValues(network, row, diffExp, gene, label+" pValue");
			updateValues(network, row, diffExp, gene, label+" FDR");
		}

		int rank = 1;
		for (String gene: rankedGenes) {
			CyRow row = queryMap.get(gene);
			if (row == null) continue;
			row.set(NAMESPACE, label+" Rank", rank);
			rank++;
		}
	}

	// Add our information to the network tables: source, experiment accession, category, category row
	public static void addNetworkColumns(ScNVManager manager, CyNetwork network) {
		CyTable netTable = network.getDefaultNetworkTable();
		createColumnIfNeeded(netTable, NAMESPACE, EXPERIMENT_SOURCE, String.class);
		createColumnIfNeeded(netTable, NAMESPACE, EXPERIMENT_ACCESSION, String.class);
		createColumnIfNeeded(netTable, NAMESPACE, CATEGORY_SOURCE, String.class);
		createColumnIfNeeded(netTable, NAMESPACE, CATEGORY, String.class);
		createColumnIfNeeded(netTable, NAMESPACE, CATEGORY_ROW, String.class);
		createListColumnIfNeeded(netTable, NAMESPACE, CATEGORY_NAMES, String.class);
	}

	public static void updateNetworkData(ScNVManager manager, CyNetwork network, 
	                                     Experiment experiment, Category category, String categoryRow) {
		CyRow row = network.getRow(network);
		
		row.set(NAMESPACE, EXPERIMENT_SOURCE, experiment.getSource().toString());
		row.set(NAMESPACE, EXPERIMENT_ACCESSION, experiment.getMetadata().get(Metadata.ACCESSION).toString());
		row.set(NAMESPACE, CATEGORY_SOURCE, category.getSource().toString());
		row.set(NAMESPACE, CATEGORY, category.toString());
		if (categoryRow != null)
			row.set(NAMESPACE, CATEGORY_ROW, categoryRow);

		List<String> categoryNames = new ArrayList<>();
		for (Object cat: experiment.getDiffExp().getLogGERMap().keySet()) {
			categoryNames.add(category.mkLabel(cat));
		}
		row.set(NAMESPACE, CATEGORY_NAMES, categoryNames);
	}

	public static Experiment getExperimentFromNetwork(ScNVManager manager, CyNetwork network) {
		if (network.getRow(network) == null || 
		    network.getDefaultNetworkTable().getColumn(NAMESPACE, EXPERIMENT_ACCESSION) == null)
			return null;

		String accession = network.getRow(network).get(NAMESPACE, EXPERIMENT_ACCESSION, String.class);
		Experiment exp = manager.getExperiment(accession);
		return exp;
	}

	public static Category getCategoryFromNetwork(ScNVManager manager, CyNetwork network) {
		Experiment exp = getExperimentFromNetwork(manager, network);
		String catName = network.getRow(network).get(NAMESPACE, CATEGORY, String.class);
		return exp.getCategory(catName);
	}

	public static String getCategoryRowFromNetwork(ScNVManager manager, CyNetwork network) {
		return network.getRow(network).get(NAMESPACE, CATEGORY_ROW, String.class);
	}

	public static List<String> getCategoryNamesFromNetwork(CyNetwork network) {
		return network.getRow(network).getList(NAMESPACE, CATEGORY_NAMES, String.class);
	}

	public static List<String> getGeneNamesFromNetwork(CyNetwork network) {
		return getGeneNamesFromNetwork(network, network.getNodeList());
	}

	public static List<String> getGeneNamesFromNetwork(CyNetwork network, List<CyNode> nodeList) {
		List<String> names = new ArrayList<>();
		for (CyNode node: nodeList) {
			names.add(getGeneNameFromNode(network, node));
		}
		return names;
	}

	public static String getGeneNameFromNode(CyNetwork network, CyNode node) {
		//TODO: add STRING namespace
		return network.getRow(node).get("query term", String.class);
	}

	public static double[] getDataFromNetwork(CyNetwork network, String column, List<CyNode> nodes) {
		double[] data = new double[nodes.size()];
		int index = 0;
		for (CyNode node: nodes) {
			Double v = network.getRow(node).get(NAMESPACE, column, Double.class);
			if (v == null)
				data[index] = Double.NaN;
			else
				data[index] = v.doubleValue();
			index++;
		}
		return data;
	}

	public static int getRowFromNode(Experiment exp, CyNetwork network, List<CyNode> nodes) {
		if (nodes == null || nodes.size() == 0 || network == null) return -1;

		// We want to use the diff exp because it's already removed the spike-ins
		List<String> rowLabels = exp.getDiffExp().getRowLabels();
		String name = getGeneNameFromNode(network, nodes.get(0));
		return rowLabels.indexOf(name);
	}

	public static List<CyNode> selectNodes(ScNVManager manager, CyNetwork network, String enrichmentType) {
		// Get DE column for this network
		Category category = getCategoryFromNetwork(manager, network);

		// First, get the current list of selected nodes
		List<CyNode> selectedNodes = CyTableUtil.getNodesInState(network, CyNetwork.SELECTED, true);

		// Now clear the selected nodes
		for (CyNode node: selectedNodes) {
			network.getRow(node).set(CyNetwork.SELECTED, false);
		}

		selectedNodes = new ArrayList<>();

		// Our row is encoded in the network name
		String categoryRow = network.getRow(network).get(CyNetwork.NAME, String.class);
		String column = categoryRow+" log2FC";
		for (CyNode node: network.getNodeList()) {
			CyRow row = network.getRow(node);
			if (row == null) continue;
			Double v = row.get(NAMESPACE, column, Double.class);
			if (v == null) {
				continue;
			}

			switch (enrichmentType) {
				case POSITIVE_ONLY:
					if (v > 0.0) {
						row.set(CyNetwork.SELECTED, true);
						selectedNodes.add(node);
					}
					break;
				case NEGATIVE_ONLY:
					if (v < 0.0) {
						row.set(CyNetwork.SELECTED, true);
						selectedNodes.add(node);
					}
					break;
			}
		}
		return selectedNodes;
	}

	/**
	 * This version of selectNodes is specifically to be used when genes are
	 * selected in a table or on a chart
	 */
	public static void selectNodes(ScNVManager manager, String accession, List<String> geneList) {
		List<CyNetwork> networks = getNetworksForAccession(manager, accession);
		Map<String, Object> args = new HashMap<>();
		for (CyNetwork net: networks) {
			args.clear();
			args.put("network", "SUID:"+net.getSUID());
			args.put("nodeList", listToString("query term:", geneList));
			manager.executeCommand("network", "select", args, true);
		}
	}

	public static void selectNodes(CyNetwork network, List<CyNode> selectedNodes) {
		for (CyNode node: selectedNodes) {
			network.getRow(node).set(CyNetwork.SELECTED, true);
		}
	}

	public static List<CyNetwork> getNetworksForAccession(ScNVManager manager, String accession) {
		CyNetworkManager netManager = manager.getService(CyNetworkManager.class);
		List<CyNetwork> nets = new ArrayList<>();

		for (CyNetwork net: netManager.getNetworkSet()) {
			if (net.getDefaultNetworkTable().getColumn(NAMESPACE, EXPERIMENT_ACCESSION) != null) {
				String acc = net.getRow(net).get(NAMESPACE, EXPERIMENT_ACCESSION, String.class);
				if (acc.equals(accession))
					nets.add(net);

			}
		}
		return nets;
	}

	public static String listToString(String prefix, List<String> list) {
		StringBuilder str = new StringBuilder();
		for (String s:list) {
			str.append(prefix+s+",");
		}
		return str.substring(0, str.length()-1);
	}

	@SuppressWarnings("unchecked")
	public static void addStyle(ScNVManager manager, CyNetwork network, String name, VisualStyle baseStyle) {
		String col = name+" log2FC";
		CyNetworkViewManager viewManager = manager.getService(CyNetworkViewManager.class);
		VisualMappingManager vizMapManager = manager.getService(VisualMappingManager.class);
		CyNetworkView view = null;
		for (CyNetworkView v: viewManager.getNetworkViews(network)) {
			if (v.getRendererId().equalsIgnoreCase("org.cytoscape.ding")) {
				view = v;
				break;
			}
		}
		if (view == null) return;

		VisualStyle newStyle = manager.getService(VisualStyleFactory.class).createVisualStyle(baseStyle);
		newStyle.setTitle(baseStyle.getTitle()+" "+name);
		newStyle.removeVisualMappingFunction(BasicVisualLexicon.NODE_FILL_COLOR);

		// Now create the new mapping
		PaletteProviderManager paletteManager = manager.getService(PaletteProviderManager.class);
		Palette rdbl = paletteManager.getPaletteProvider("ColorBrewer").getPalette("Red-Blue", 9);
		Color[] colors = rdbl.getColors();

		// Get a function factory
		VisualMappingFunctionFactory vmff = manager.getService(VisualMappingFunctionFactory.class, 
		                                                       "(mapping.type=continuous)");

		@SuppressWarnings("rawtypes")
		ContinuousMapping colorMapping =
		   (ContinuousMapping) vmff.createVisualMappingFunction("scNetViz::"+col, 
		                                                        Double.class, BasicVisualLexicon.NODE_FILL_COLOR);

		double[] minMax = findMinMax(network, col);
		colorMapping.addPoint (minMax[0], new BoundaryRangeValues<Paint>(colors[8], colors[7], colors[7]));
		colorMapping.addPoint (0, new BoundaryRangeValues<Paint>(colors[4], colors[4], colors[4]));
		colorMapping.addPoint (minMax[1], new BoundaryRangeValues<Paint>(colors[1], colors[1], colors[0]));

		newStyle.addVisualMappingFunction(colorMapping);
		vizMapManager.setVisualStyle(newStyle, view);
		// vizMapManager.addVisualStyle(newStyle);
	}

	public static VisualStyle getVisualStyle(ScNVManager manager, String name) {
		VisualMappingManager vizMapManager = manager.getService(VisualMappingManager.class);
		for (VisualStyle style: vizMapManager.getAllVisualStyles()) {
			if (style.getTitle().startsWith(name))
				return style;
		}
		return null;
	}

	public static Map<String, CyRow> createQueryMap(CyNetwork network) {
		Map<String, CyRow> qMap = new HashMap<>();
		for (CyNode node: network.getNodeList()) {
			qMap.put(network.getRow(node).get("query term", String.class), network.getRow(node));
		}
		return qMap;
	}

	public static List<CyNode> getSelectedNodes(CyNetwork network) {
		return CyTableUtil.getNodesInState(network, CyNetwork.SELECTED, true);
	}

	public static double[] findMinMax(CyNetwork network, String col) {
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		// System.out.println("Looking at: "+col);
		for (CyNode node: network.getNodeList()) {
			Double v = network.getRow(node).get(NAMESPACE, col, Double.class);
			if (v != null && !Double.isNaN(v)) {
				if (v < min) min = v;
				if (v > max) max = v;
			}
		}
		if (Math.abs(min) > max)
			max = Math.abs(min);
		else if (min < 0)
			min = -max;
		else
			min = 0;

		double[] res = new double[2];
		res[0] = min;
		res[1] = max;
		// System.out.println("min = "+min+", max = "+max);
		return res;
	}

	public static void updateValues(CyNetwork network, CyRow row, DoubleMatrix dMat, 
	                                String rowLabel, String colLabel) {
		double v = dMat.getDoubleValue(rowLabel, colLabel);
		if (!Double.isNaN(v) && !Double.isInfinite(v))
			row.set(NAMESPACE, colLabel, v);
	}

	public static void createColumnIfNeeded(CyTable table, String namespace, String column, Class<?> clazz) {
		if (table.getColumn(namespace, column) == null)
			table.createColumn(namespace, column, clazz, false);
	}

	public static void createListColumnIfNeeded(CyTable table, String namespace, String column, Class<?> clazz) {
		if (table.getColumn(namespace, column) == null)
			table.createListColumn(namespace, column, clazz, false);
	}

	public static boolean isInteger(String v) {
		try {
			int a = Integer.parseInt(v.trim());
			return true;
		} catch (NumberFormatException e) {
			return false;
		}
	}
}
