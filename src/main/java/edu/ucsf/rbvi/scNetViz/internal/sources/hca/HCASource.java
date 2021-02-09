package edu.ucsf.rbvi.scNetViz.internal.sources.hca;

import java.io.File;
import java.net.URLEncoder;
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
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.tasks.HCAListEntriesTaskFactory;
// import edu.ucsf.rbvi.scNetViz.internal.sources.hca.tasks.HCALoadExperimentTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.tasks.HCAShowEntriesTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.view.HCAEntryFrame;

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

public class HCASource implements Source {
	public static String HCA_PROJECT_URL = "https://service.azul.data.humancellatlas.org/index/projects?catalog=dcp1&";

	final Logger logger;
	final ScNVManager scNVManager;
	final Map<String, HCAMetadata> metadataMap;
	HCAEntryFrame entryFrame = null;

	public HCASource (ScNVManager manager) {
		scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);

		// Read in all of the metadata (entries) for hca
		metadataMap = new HashMap<>();

		{
			Properties props = new Properties();
			props.put(SOURCENAME, "hca");
			scNVManager.registerService(this, Source.class, props);
		}

		// Register our task factories
		{
			Properties props = new Properties();
			props.setProperty(TITLE, "From the Human Cell Atlas...");
			props.setProperty(PREFERRED_MENU, "Apps.scNetViz.Load Experiment[10.1]");
			props.setProperty(IN_TOOL_BAR, "TRUE");
			props.setProperty(TOOL_BAR_GRAVITY, "100.0");
			props.setProperty(MENU_GRAVITY, "20.0");
			props.setProperty(TOOLTIP, "Browse the Human Cell Atlas");
			String ebiLogoURL = getClass().getResource("/images/hca.png").toString();
			props.setProperty(LARGE_ICON_URL, ebiLogoURL);
			scNVManager.registerService(new HCAShowEntriesTaskFactory(manager, this), TaskFactory.class, props);
		}

		// Register our commands
		{
			Properties props = new Properties();
			props.setProperty(COMMAND_DESCRIPTION, "List all Gene Expression Atlas (HCA) entries available");
			props.setProperty(COMMAND_NAMESPACE, "scnetviz");
			props.setProperty(COMMAND, "list hca entries");
			props.setProperty(COMMAND_SUPPORTS_JSON, "true");
			props.setProperty(COMMAND_EXAMPLE_JSON, "{}");

			scNVManager.registerService(new HCAListEntriesTaskFactory(manager, this), TaskFactory.class, props);
		}
		
		{
			Properties props = new Properties();
			props.put(COMMAND_NAMESPACE, "scnetviz");
			props.put(COMMAND_DESCRIPTION, "Load an experiment from the Human Cell Atlas matrix service");
			props.put(COMMAND, "load hca experiment");
			props.put(COMMAND_SUPPORTS_JSON, "true");
			// scNVManager.registerService(new HCALoadExperimentTaskFactory(manager, this), TaskFactory.class, props);
		}

	}

	public String getName() { return "HCA"; }
	public String toString() { return "Human Cell Atlas"; }

	public void loadHCAEntries(TaskMonitor taskMonitor) {
		if (metadataMap.size() > 0) return;
		// String query = "filters="+URLEncoder.encode("{\"file\":{\"fileFormat\":{\"is\":[\"matrix\"]}}}");
		String query = "filters="+URLEncoder.encode("{\"fileFormat\":{\"is\":[\"matrix\"]}}");
		try {
			JSONObject json = HTTPUtils.fetchJSON(HCA_PROJECT_URL+query, taskMonitor);
			JSONArray hits = (JSONArray) json.get("hits");
      while (json.containsKey("pagination") && 
             ((JSONObject)json.get("pagination")).containsKey("next") && 
             ((JSONObject)json.get("pagination")).get("next") != null) {
			  json = HTTPUtils.fetchJSON(((JSONObject)json.get("pagination")).get("next").toString(), taskMonitor);
        hits.addAll((JSONArray)json.get("hits"));
      }
			for (Object hit: hits) {
				// For some reason, even though we only ask for projects with matrix files, we get others also
				if (HCAMetadata.hasMatrix((JSONObject)hit)) {
          List<String> organs = HCAMetadata.getOrgans((JSONObject)hit);
          for (String organ: organs) {
					  HCAMetadata entry = new HCAMetadata(organ, (JSONObject)hit);
            if (entry.get(HCAMetadata.MATRIX) != null)
					    metadataMap.put((String)entry.get(Metadata.ACCESSION), entry);
          }
				}
			}
		} catch (Exception e) {
      e.printStackTrace();
			taskMonitor.showMessage(TaskMonitor.Level.ERROR, "Exception reading HCA experiment list: "+e.getMessage());
		}
  }

	public List<String> getAccessions() {
		if (metadataMap.size() == 0)
			loadHCAEntries(null);
		return new ArrayList<>(metadataMap.keySet());
	}

	public List<Metadata> getMetadata() {
		if (metadataMap.size() == 0)
			loadHCAEntries(null);
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
		HCAExperiment exp = new HCAExperiment(scNVManager, this, (HCAMetadata)metadata);
		exp.fetchMTX (monitor);
		exp.fetchDesign(monitor);
		return exp;
	}

	public HCAEntryFrame getEntryFrame() { return entryFrame; }

	public void showEntriesTable(boolean showHide) {
		if (entryFrame == null && showHide) {
			entryFrame = new HCAEntryFrame(scNVManager, this);
		} else if (showHide) {
			entryFrame.setVisible(true);
		} else if (entryFrame != null)
			entryFrame.setVisible(false);
		
	}

	public Experiment loadExperimentFromSession(JSONObject jsonExperiment, Map<String, File> fileMap) {
		return null;
	}

	public Category loadCategoryFromSession(JSONObject jsonExperiment, Experiment experiment, Map<String, File> fileMap) {
		return null;
	}

	public DifferentialExpression loadDiffExpFromSession(JSONObject jsonDiffExp, Experiment experiment, Map<String, File> fileMap) {
		return null;
	}
}
