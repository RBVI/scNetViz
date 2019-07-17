package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import org.apache.log4j.Logger;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks.GXAListEntriesTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks.GXALoadExperimentTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks.GXAShowEntriesTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.view.GXAEntryFrame;

import static org.cytoscape.work.ServiceProperties.COMMAND;
import static org.cytoscape.work.ServiceProperties.COMMAND_DESCRIPTION;
import static org.cytoscape.work.ServiceProperties.COMMAND_NAMESPACE;
import static org.cytoscape.work.ServiceProperties.COMMAND_SUPPORTS_JSON;
import static org.cytoscape.work.ServiceProperties.COMMAND_EXAMPLE_JSON;

import static org.cytoscape.work.ServiceProperties.ID;
import static org.cytoscape.work.ServiceProperties.IN_MENU_BAR;
import static org.cytoscape.work.ServiceProperties.IN_TOOL_BAR;
import static org.cytoscape.work.ServiceProperties.INSERT_SEPARATOR_BEFORE;
import static org.cytoscape.work.ServiceProperties.LARGE_ICON_URL;
import static org.cytoscape.work.ServiceProperties.MENU_GRAVITY;
import static org.cytoscape.work.ServiceProperties.PREFERRED_MENU;
import static org.cytoscape.work.ServiceProperties.TITLE;
import static org.cytoscape.work.ServiceProperties.TOOL_BAR_GRAVITY;
import static org.cytoscape.work.ServiceProperties.TOOLTIP;

public class GXASource implements Source {
	public static String EXPERIMENTS_URL = "https://www.ebi.ac.uk/gxa/sc/json/experiments";

	final Logger logger;
	final ScNVManager scNVManager;
	final Map<String, GXAMetadata> metadataMap;
	GXAEntryFrame entryFrame = null;

	public GXASource (ScNVManager manager) {
		scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);

		// Read in all of the metadata (entries) for gxa
		metadataMap = new HashMap<>();

		{
			Properties props = new Properties();
			props.put(SOURCENAME, "gxa");
			scNVManager.registerService(this, Source.class, props);
		}

		// Register our task factories
		{
			Properties props = new Properties();
			props.put(TITLE, "Browse Single Cell Expression Atlas");
			props.put(PREFERRED_MENU, "Apps.scNetViz");
			props.setProperty(IN_TOOL_BAR, "TRUE");
			props.setProperty(TOOL_BAR_GRAVITY, "100f");
			props.setProperty(TOOLTIP, "Show Experiments Table");
			String ebiLogoURL = getClass().getResource("/images/EMBL-EBI-Logo-36x36.png").toString();
			props.setProperty(LARGE_ICON_URL, ebiLogoURL);
			scNVManager.registerService(new GXAShowEntriesTaskFactory(manager, this), TaskFactory.class, props);
		}

		// Register our commands
		{
			Properties props = new Properties();
			props.setProperty(COMMAND_DESCRIPTION, "List all Gene Expression Atlas (GXA) entries available");
			props.setProperty(COMMAND_NAMESPACE, "scnetviz");
			props.setProperty(COMMAND, "list gxa entries");
			props.setProperty(COMMAND_SUPPORTS_JSON, "true");
			props.setProperty(COMMAND_EXAMPLE_JSON, "{}");

			scNVManager.registerService(new GXAListEntriesTaskFactory(manager, this), TaskFactory.class, props);
		}
		
		{
			Properties props = new Properties();
			props.put(COMMAND_NAMESPACE, "scnetviz");
			props.put(COMMAND_DESCRIPTION, "Load an experiment from the EBI Gene Expression Atlas");
			props.put(COMMAND, "load gxa experiment");
			props.put(COMMAND_SUPPORTS_JSON, "true");
			scNVManager.registerService(new GXALoadExperimentTaskFactory(manager, this), TaskFactory.class, props);
		}

	}

	public String getName() { return "GXA"; }
	public String toString() { return "Single Cell Expression Atlas"; }

	public void loadGXAEntries(TaskMonitor taskMonitor) {
		if (metadataMap.size() > 0) return;
		JSONObject json;
		try {
			json = HTTPUtils.fetchJSON(EXPERIMENTS_URL, taskMonitor);
		} catch (Exception e) {
			taskMonitor.showMessage(TaskMonitor.Level.ERROR, "Unable to fetch expreiments from EBI: "+
			                        e.getMessage());
			return;
		}
		if (json == null) return;
		JSONArray experiments = (JSONArray) json.get("aaData");
		for (Object exp: experiments) {
			GXAMetadata entry = new GXAMetadata((JSONObject)exp);
			metadataMap.put((String)entry.get(Metadata.ACCESSION), entry);
		}
  }

	public List<String> getAccessions() {
		if (metadataMap.size() == 0)
			loadGXAEntries(null);
		return new ArrayList<>(metadataMap.keySet());
	}

	public List<Metadata> getMetadata() {
		if (metadataMap.size() == 0)
			loadGXAEntries(null);
		return new ArrayList<>(metadataMap.values());
	}

	@Override
	public Experiment getExperiment(String accession) {
		return getExperiment(metadataMap.get(accession), null, false);
	}

	public Experiment getExperiment(String accession, boolean showTable) {
		return getExperiment(metadataMap.get(accession), null, showTable);
	}

	public Experiment getExperiment(Metadata metadata, boolean showTable) {
		return getExperiment(metadata, null, showTable);
	}

	public Experiment getExperiment(String accession, TaskMonitor monitor, boolean showTable) {
		return getExperiment(metadataMap.get(accession), monitor, showTable);
	}

	public Experiment getExperiment(Metadata metadata, TaskMonitor monitor, boolean showTable) {
		GXAExperiment exp = new GXAExperiment(scNVManager, this, (GXAMetadata)metadata);
		exp.fetchMTX (monitor);
		if (showTable) {
			exp.fetchClusters(monitor);
			exp.fetchDesign(monitor);
		} else {
			exp.fetchClusters();
			exp.fetchDesign();
		}
		return exp;
	}

	public GXAEntryFrame getEntryFrame() { return entryFrame; }

	public void showEntriesTable(boolean showHide) {
		if (entryFrame == null && showHide) {
			entryFrame = new GXAEntryFrame(scNVManager, this);
		} else if (showHide) {
			entryFrame.setVisible(true);
		} else if (entryFrame != null)
			entryFrame.setVisible(false);
	}

	public Experiment loadExperimentFromSession(JSONObject jsonExperiment, Map<String, File> fileMap) {
		// Get and create the GXAMetadata from the JSON
		GXAMetadata metadata = new GXAMetadata();
		if (jsonExperiment.containsKey("metadata"))
			metadata.fromJSON((JSONObject)jsonExperiment.get("metadata"));

		// Create the GXAExperiment
		GXAExperiment experiment = new GXAExperiment(scNVManager, this, metadata);
		try {
			experiment.loadFromSession(fileMap);
		} catch (Exception e) {
			logger.error("Unable to load experiment from session: "+e.toString());
			return null;
		}

		return experiment;
	}

	public Category loadCategoryFromSession(JSONObject jsonCategory, Experiment experiment, Map<String, File> fileMap) {
		if (experiment instanceof GXAExperiment) {
			Category category = null;
			try {
				category = ((GXAExperiment)experiment).loadCategoryFromSession(jsonCategory, fileMap);
			} catch (Exception e) {
				logger.error("Unable to load category from session: "+e.toString());
				return null;
			}
			return category;
		}
		return null;
	}
	public DifferentialExpression loadDiffExpFromSession(JSONObject jsonDiffExp, Experiment experiment, Map<String, File> fileMap) {
		try {
			return ((GXAExperiment)experiment).loadDiffExpFromSession(jsonDiffExp, fileMap);
		} catch (Exception e) {
			logger.error("Unable to load differential expression from session: "+e.toString());
			return null;
		}

	}
}
