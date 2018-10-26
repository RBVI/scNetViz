package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;
import org.apache.log4j.Logger;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;

public class GXASource implements Source {
	public static String EXPERIMENTS_URL = "https://www.ebi.ac.uk/gxa/sc/json/experiments";

	final Logger logger;
	final ScNVManager scNVManager;
	final Map<String, GXAMetadata> metadataMap;

	public GXASource (ScNVManager manager) {
		scNVManager = manager;
		logger = Logger.getLogger(CyUserLog.NAME);

		// Read in all of the metadata (entries) for gxa
		metadataMap = new HashMap<>();
	}

	public String getSourceName() { return "EBI GXA"; }

	public void loadGXAEntries(TaskMonitor taskMonitor) {
		if (metadataMap.size() > 0) return;
		JSONObject json = HTTPUtils.fetchJSON(EXPERIMENTS_URL, taskMonitor);
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

	public Experiment getExperiment(String accession) {
		return getExperiment(accession, null);
	}

	public Experiment getExperiment(Metadata metadata) {
		return getExperiment(metadata, null);
	}

	public Experiment getExperiment(Metadata metadata, TaskMonitor monitor) {
		return getExperiment((String)metadata.get(Metadata.ACCESSION), monitor);
	}

	public Experiment getExperiment(String accession, TaskMonitor monitor) {
		GXAExperiment exp = new GXAExperiment(scNVManager, this);
		exp.fetchMTX (accession, monitor);
		exp.fetchClusters();
		exp.fetchDesign();
		return exp;
	}
}
