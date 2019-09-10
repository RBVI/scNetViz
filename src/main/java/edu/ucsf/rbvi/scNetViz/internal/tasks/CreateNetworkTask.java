package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.event.CyEventHelper;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.layout.CyLayoutAlgorithmManager;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.View;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.TaskObserver;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.TunableSetter;
import org.cytoscape.work.json.JSONResult;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

// Tunable to choose experiment?

public class CreateNetworkTask extends AbstractTask implements ObservableTask {
	final ScNVManager manager;
	final CyEventHelper cyEventHelper;
	CyNetwork unionNetwork = null;
	VisualStyle baseStyle = null;
	Experiment experiment;
	final CyApplicationManager appManager;

	final DifferentialExpression diffExp;
	double fdrCutoff;
	double log2FCCutoff;
	int topGenes;
	int maxGenes;
	boolean positiveOnly;

	public CreateNetworkTask(final ScNVManager manager, DifferentialExpression diffExp, 
	                         double fdrCutoff, double log2FCCutoff, int nGenes, boolean positiveOnly,
	                         int maxGenes) {
		super();
		this.manager = manager;
		this.diffExp = diffExp;
		this.fdrCutoff = fdrCutoff;
		this.log2FCCutoff = log2FCCutoff;
		this.topGenes = nGenes;
		this.maxGenes = maxGenes;
		this.positiveOnly = positiveOnly;
		cyEventHelper = manager.getService(CyEventHelper.class);
		appManager = manager.getService(CyApplicationManager.class);
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("Creating Networks");

		Category category = diffExp.getCurrentCategory();
		experiment = category.getExperiment();
		Set<Object> categoryValues = diffExp.getCategoryValues();
		List<String> rowLabels = category.getMatrix().getRowLabels();
		String categoryRow = category.toString()+" ("+rowLabels.get(category.getSelectedRow())+")";

		Map<Object, List<String>> geneMap = new HashMap<>();

		// Use a hash set to avoid duplicate genes
		Set<String> allGenes = new HashSet<String>();

		// Iterate over each category value
		for (Object cat: categoryValues) {
			// Get the genes that match our criteria
			List<String> geneList = new ArrayList<>();
			diffExp.getGeneList(cat, fdrCutoff, log2FCCutoff, 
			                    topGenes, positiveOnly, maxGenes, geneList);
			if (geneList != null && geneList.size() > 0) {
				allGenes.addAll(geneList);
				geneMap.put(cat, geneList);
				// System.out.println("Ordered gene list for "+cat+" "+geneList);
			} else {
				monitor.showMessage(TaskMonitor.Level.WARN, "No genes passed the cutoff for "+category.mkLabel(cat));
			}
		}

		if (allGenes.size() == 0) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "No genes passed the cutoffs for any cluster");
			return;
		}

		// Create the union network
		// Create the network
		createStringNetwork(null, categoryRow, new ArrayList<String>(allGenes), geneMap, monitor);

		for (Object cat: categoryValues) {
			List<String> geneList = geneMap.get(cat);
			createSubNetwork(cat, category.mkLabel(cat), geneList, monitor);
		}
	}

	// TODO: return the networks?
	public <R> R getResults(Class<? extends R> clazz) {
		if (clazz.equals(DifferentialExpression.class))
			return (R)diffExp;
		return null;
	}

	private void createStringNetwork(Object cat, String name, List<String> geneList, 
	                                 Map<Object, List<String>> geneMap, TaskMonitor monitor) {
		monitor.setTitle("Retrieving STRING network for: "+name);
		monitor.showMessage(TaskMonitor.Level.INFO, "Retrieving STRING network for: "+name);
		Map<String, Object> args = new HashMap<>();
		args.put("query", listToString(geneList, ""));
		args.put("species", diffExp.getCurrentCategory().getExperiment().getSpecies());
		args.put("limit", "0");
		args.put("includesViruses", "false");
		manager.executeCommand("string", "protein query", args, 
		                       new RenameNetwork(diffExp, cat, name, geneList, geneMap, false, monitor), true);
		cyEventHelper.flushPayloadEvents();
	}

	private void createSubNetwork(Object cat, String name, List<String> geneList, TaskMonitor monitor) {
		monitor.setTitle("Creating subnetwork for: "+name);
		monitor.showMessage(TaskMonitor.Level.INFO, "Creating subnetwork for: "+name);
		Map<String, Object> args = new HashMap<>();
		args.put("nodeList", listToString(geneList, "query term:"));
		args.put("networkName", name);
		args.put("source", "SUID:"+unionNetwork.getSUID());
		manager.executeCommand("network", "create", args,
		                       new RenameNetwork(diffExp, cat, name, geneList, null, false, monitor), true);
	}

	private String listToString(List<String> list, String prefix) {
		if (list == null || list.size() < 1) return "";
		String str = prefix+list.get(0);
		for (int i = 1; i < list.size(); i++) {
			str += ","+prefix+list.get(i);
		}
		return str;
	}

	private class RenameNetwork implements TaskObserver {
		String name;
		List<String> geneList;
		final DifferentialExpression diffExp;
		boolean getEnrichment = false;
		Object cat = null;
		final TaskMonitor monitor;
		final Map<Object, List<String>> geneMap;

		public RenameNetwork(final DifferentialExpression diffExp, Object cat, String newName, 
		                     final List<String> geneList, final Map<Object, List<String>> geneMap,
		                     boolean getEnrichment, final TaskMonitor monitor) {
			this.name = newName;
			this.cat = cat;
			this.geneList = geneList;
			this.diffExp = diffExp;
			this.monitor = monitor;
			this.getEnrichment = getEnrichment;
			this.geneMap = geneMap;
		}

		public void allFinished(FinishStatus status) {}
		
		public void taskFinished(ObservableTask task) {
			Object res = task.getResults(JSONResult.class);
			if (res == null) return;
			if (task.toString().contains("RegisterNetworkTask"))
				return;

			// Get some information we're going to need later
			Category category = diffExp.getCurrentCategory();
			Experiment experiment = category.getExperiment();
			List<String> rowLabels = category.getMatrix().getRowLabels();
			String categoryRow = rowLabels.get(category.getSelectedRow());

			CyNetwork network = null;
			CyNetworkView networkView = null;
			if (res instanceof JSONResult) {
				network = ModelUtils.getNetworkFromJSON(manager, (JSONResult)res);
				// System.out.println("Got network "+network+" from "+((JSONResult)res).getJSON());
			} else if (res instanceof CyNetworkView) {
				networkView = (CyNetworkView)res;
				network = ((CyNetworkView)res).getModel();
			}
			if (cat == null) {
				ModelUtils.rename(network, network, name+" Network");
				cyEventHelper.flushPayloadEvents();
			}

			baseStyle = ModelUtils.getVisualStyle(manager, "STRING style");

			Map<String, Object> args = new HashMap<>();
			args.put("network", network.getRow(network).get(CyNetwork.NAME, String.class));
			manager.executeCommand("view", "set current", args, null, true);
			cyEventHelper.flushPayloadEvents();

			manager.executeCommand("string", "make string", args, null, true);

			manager.executeCommand("string", "hide images", args, null, true);
			manager.executeCommand("string", "hide glass", args, null, true);

			if (getEnrichment) {
				monitor.setTitle("Retrieving enrichment for : "+name);
				args.clear();

				manager.executeCommand("string", "retrieve enrichment", args, null, true);
				manager.executeCommand("string", "show enrichment", args, null, true);
				manager.executeCommand("string", "show charts", args, null, true);
			}

			// Create the columns
			monitor.setTitle("Adding data to network for: "+name);
			if (cat != null) {
				// Style the network
				monitor.setTitle("Creating network style for: "+name);
				ModelUtils.addStyle(manager, network,  name, baseStyle);
				// Add our information to the network tables: source, experiment accession, category, category row
				ModelUtils.addNetworkColumns(manager, network);
				ModelUtils.updateNetworkData(manager, network, experiment, category, categoryRow);

				// And re-lay it out
				CyLayoutAlgorithm alg = manager.getService(CyLayoutAlgorithmManager.class).getLayout("force-directed");
				Object context = alg.createLayoutContext();
				TunableSetter setter = manager.getService(TunableSetter.class);
				Map<String, Object> layoutArgs = new HashMap<>();
				layoutArgs.put("defaultNodeMass", 10.0);
				setter.applyTunables(context, layoutArgs);
				Set<View<CyNode>> nodeViews = new HashSet<>(networkView.getNodeViews());

				Task resultsPanelTask = new ShowResultsPanelTask(manager, experiment);
				TaskIterator ti = new TaskIterator();
				ti.append(alg.createTaskIterator(networkView, context, nodeViews, "score"));
				ti.append(resultsPanelTask);
				manager.executeTasks(ti);
			} else {
				monitor.setTitle("Adding data to network for: "+name);
				unionNetwork = network;
				Set<Object> categoryValues = diffExp.getCategoryValues();
				// Add our information to the network tables: source, experiment accession, category
				ModelUtils.addNetworkColumns(manager, network);
				ModelUtils.updateNetworkData(manager, network, experiment, category, null);
				for (Object cat1: categoryValues) {
					List<String> sortedList = geneMap.get(cat1);

					ModelUtils.createDEColumns(manager, network, diffExp, category.mkLabel(cat1));
					// Add the data
					ModelUtils.updateDEData(manager, network, geneList, diffExp, category.mkLabel(cat1), sortedList);

				}
				// Task resultsPanelTask = new ShowResultsPanelTask(manager, experiment);
				// manager.executeTasks(new TaskIterator(resultsPanelTask));
			}
		}
	}
}
