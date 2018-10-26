package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;

public class GXAMetadata extends HashMap<String, Object> implements Metadata {
	public static String KINGDOM = "kingdom";
	public static String CONTRASTS = "contrasts";
	public static String ASSAYS = "assays";
	public static String FACTORS = "factors";

	public GXAMetadata(JSONObject json) {
		super();
		put(TYPE,(String) json.get("experimentType"));
		put(ACCESSION,(String) json.get("experimentAccession"));
		put(DESCRIPTION, (String) json.get("experimentDescription"));
		put(DATE, (String) json.get("lastUpdate"));
		put(ASSAYS, (Long) json.get("numberOfAssays"));
		put(CONTRASTS, (Long) json.get("numberOfContrasts"));
		put(SPECIES, (String) json.get("species"));
		put(KINGDOM, (String) json.get("kingdom"));
		JSONArray factors = (JSONArray) json.get("experimentalFactors");
		List<String> expFactors = new ArrayList<String>();
		for (Object obj: factors) {
			expFactors.add((String)obj);
		}
		put(FACTORS, expFactors);
	}

	public String toString() {
		return "<html><p style='width: 500px'><b>"+get(ACCESSION)+"</b>: "+get(DESCRIPTION)+"</p></html>";
	}
}
