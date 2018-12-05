package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.awt.Color;
import java.awt.Paint;
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
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class ModelUtils {
	public static CyNetwork getNetworkFromJSON(ScNVManager manager, JSONResult res) {
		try {
			JSONObject json = (JSONObject) new JSONParser().parse(res.getJSON());
			long networkSUID = ((Number)json.get("SUID")).longValue();
			CyNetwork network = manager.getService(CyNetworkManager.class).getNetwork(networkSUID);
			return network;
		} catch (Exception e) { return null; }
	}

	public static void rename(CyNetwork network, CyIdentifiable id, String name) {
		network.getRow(id).set(CyNetwork.NAME, name);
	}

	public static void createDEColumns(ScNVManager manager, CyNetwork network, DifferentialExpression diffExp, 
	                                   String label) {
		CyTable nodeTable = network.getDefaultNodeTable();
		nodeTable.createColumn("scNetViz", label+" MTC", Double.class, false);
		nodeTable.createColumn("scNetViz", label+" Min.pct", Double.class, false);
		nodeTable.createColumn("scNetViz", label+" MDTC", Double.class, false);
		nodeTable.createColumn("scNetViz", label+" logGER", Double.class, false);
		nodeTable.createColumn("scNetViz", label+" pValue", Double.class, false);
		nodeTable.createColumn("scNetViz", label+" FDR", Double.class, false);
	}

	public static void updateDEData(ScNVManager manager, CyNetwork network, 
	                                List<String> geneList, DifferentialExpression diffExp, 
	                                String label) {
		Map<String, CyRow> queryMap = createQueryMap(network);
		for (String gene: geneList) {
			CyRow row = queryMap.get(gene);
			if (row == null) continue;

			updateValues(network, row, diffExp, gene, label+" MTC");
			updateValues(network, row, diffExp, gene, label+" Min.pct");
			updateValues(network, row, diffExp, gene, label+" MDTC");
			updateValues(network, row, diffExp, gene, label+" logGER");
			updateValues(network, row, diffExp, gene, label+" pValue");
			updateValues(network, row, diffExp, gene, label+" FDR");
		}
	}

	@SuppressWarnings("unchecked")
	public static void addStyle(ScNVManager manager, CyNetwork network, String name) {
		String col = name+" logGER";
		System.out.println("Column = "+col);
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

		VisualStyle style = vizMapManager.getVisualStyle(view);
		VisualStyle newStyle = manager.getService(VisualStyleFactory.class).createVisualStyle(style);
		newStyle.setTitle(style.getTitle()+" "+name);
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
		   (ContinuousMapping) vmff.createVisualMappingFunction("scNetViz::"+col, Double.class, BasicVisualLexicon.NODE_FILL_COLOR);

		double[] minMax = findMinMax(network, col);
		colorMapping.addPoint (minMax[0], new BoundaryRangeValues<Paint>(colors[0], colors[1], colors[1]));
		colorMapping.addPoint (0, new BoundaryRangeValues<Paint>(colors[4], colors[4], colors[4]));
		colorMapping.addPoint (minMax[1], new BoundaryRangeValues<Paint>(colors[7], colors[7], colors[8]));

		newStyle.addVisualMappingFunction(colorMapping);
		vizMapManager.addVisualStyle(newStyle);
		vizMapManager.setVisualStyle(newStyle, view);
	}

	public static Map<String, CyRow> createQueryMap(CyNetwork network) {
		Map<String, CyRow> qMap = new HashMap<>();
		for (CyNode node: network.getNodeList()) {
			qMap.put(network.getRow(node).get("query term", String.class), network.getRow(node));
		}
		return qMap;
	}

	public static double[] findMinMax(CyNetwork network, String col) {
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		System.out.println("Looking at: "+col);
		for (CyNode node: network.getNodeList()) {
			Double v = network.getRow(node).get("scNetViz", col, Double.class);
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
		System.out.println("min = "+min+", max = "+max);
		return res;
	}

	public static void updateValues(CyNetwork network, CyRow row, DoubleMatrix dMat, 
	                                String rowLabel, String colLabel) {
		double v = dMat.getDoubleValue(rowLabel, colLabel);
		if (!Double.isNaN(v) && !Double.isInfinite(v))
			row.set("scNetViz", colLabel, v);
	}
}
