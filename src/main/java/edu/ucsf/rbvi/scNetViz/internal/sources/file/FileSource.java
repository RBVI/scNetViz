package edu.ucsf.rbvi.scNetViz.internal.sources.file;

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

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ShowExperimentTableTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileExperimentTaskFactory;

import static org.cytoscape.work.ServiceProperties.COMMAND;
import static org.cytoscape.work.ServiceProperties.COMMAND_DESCRIPTION;
import static org.cytoscape.work.ServiceProperties.COMMAND_NAMESPACE;
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

public class FileSource implements Source {
	final Logger logger;
	final ScNVManager scNVManager;
	final Map<String, FileMetadata> metadataMap;

	public FileSource (ScNVManager manager) {
		scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);

		// Read in all of the metadata (entries) for gxa
		metadataMap = new HashMap<>();

		{
			Properties props = new Properties();
			props.put(SOURCENAME, "file");
			scNVManager.registerService(this, Source.class, props);
		}

		// Register our task factories
		{
			Properties props = new Properties();
			props.put(TITLE, "Import Experiment from File");
			props.put(PREFERRED_MENU, "Apps.scNetViz");
			props.setProperty(IN_TOOL_BAR, "FALSE");
			props.setProperty(IN_MENU_BAR, "TRUE");
			props.put(COMMAND_NAMESPACE, "scnetviz");
			props.put(COMMAND_DESCRIPTION, "Load an experiment from a file");
			props.put(COMMAND, "load experiment");
			scNVManager.registerService(new FileExperimentTaskFactory(manager, this), TaskFactory.class, props);
		}

		{
			Properties props = new Properties();
			props.put(TITLE, "Add category to an experiment");
			props.put(PREFERRED_MENU, "Apps.scNetViz");
			props.setProperty(IN_TOOL_BAR, "FALSE");
			props.setProperty(IN_MENU_BAR, "TRUE");
			props.put(COMMAND_NAMESPACE, "scnetviz");
			props.put(COMMAND_DESCRIPTION, "Add a category to an experiment from a file");
			props.put(COMMAND, "load category");
			scNVManager.registerService(new FileCategoryTaskFactory(manager, this), TaskFactory.class, props);
		}
		
	}

	public String getName() { return "File"; }
	public String toString() { return "File Experiment"; }

	public List<Metadata> getMetadata() {
		return new ArrayList<>(metadataMap.values());
	}

	public List<String> getAccessions() {
		return new ArrayList<>(metadataMap.keySet());
	}

	public Experiment getExperiment(String accession) {
		return getExperiment(metadataMap.get(accession), null);
	}

	public Experiment getExperiment(FileMetadata metadata) {
		return getExperiment(metadata, null);
	}

	public Experiment getExperiment(String accession, TaskMonitor monitor) {
		return getExperiment(metadataMap.get(accession), monitor);
	}

	public Experiment getExperiment(FileMetadata metadata, TaskMonitor monitor) {
		FileExperiment exp = new FileExperiment(scNVManager, this, metadata);
		if (exp != null) {
			exp.readMTX(monitor);
			metadataMap.put((String)metadata.get(Metadata.ACCESSION), metadata);
		}
		return exp;
	}
}
