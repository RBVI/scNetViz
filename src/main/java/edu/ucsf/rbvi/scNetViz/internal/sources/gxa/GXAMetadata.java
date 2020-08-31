package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import edu.ucsf.rbvi.scNetViz.internal.utils.JSONUtils;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;

public class GXAMetadata extends HashMap<String, Object> implements Metadata {
	public static String KINGDOM = "kingdom";
	public static String CONTRASTS = "contrasts";
	public static String ASSAYS = "assays";
	public static String FACTORS = "factors";
	public static String TECHTYPE = "technologyType";
	public static String PROJECTS = "projects";

	public GXAMetadata(JSONObject json) {
		super();
		// System.out.println("Metadata: "+json.toString());
		safePut(TYPE,(String) json.get("rawExperimentType"));
		safePut(ACCESSION,(String) json.get("experimentAccession"));
		safePut(DESCRIPTION, (String) json.get("experimentDescription"));
		safePut(DATE, (String) json.get("lastUpdate"));
		safePut(ASSAYS, (Long) json.get("numberOfAssays"));
		safePut(CONTRASTS, (Long) json.get("numberOfContrasts"));
		safePut(SPECIES, (String) json.get("species"));
		safePut(KINGDOM, (String) json.get("kingdom"));
		safePut(TECHTYPE, JSONUtils.jsonArrayToList((JSONArray) json.get(TECHTYPE), String.class));
		safePut(PROJECTS, JSONUtils.jsonArrayToList((JSONArray) json.get("experimentProjects"), String.class));
		JSONArray factors = (JSONArray) json.get("experimentalFactors");
		safePut(FACTORS, JSONUtils.jsonArrayToList(factors, String.class));
	}

	public GXAMetadata() {
		super();
	}

	public String toHTML() {
		return "<html><p style='width: 500px'><b>"+get(ACCESSION)+"</b>: "+get(DESCRIPTION)+"</p></html>";
	}

	public String toString() {
		return get(ACCESSION)+": "+get(DESCRIPTION);
	}

	public String toJSON() {
		String json = "{";
		for (String key: keySet()) {
			Object v = get(key);
			if (v == null) {
				System.out.println("No entry for "+key);
				continue;
			}
			json += "\""+key+"\":";
			if (v instanceof List) {
				List vL = (List)v;
				// Check for empty list
				if (vL.isEmpty()) {
					json += "[],";
				} else {
					json += "[";
					for (Object le: vL) {
						json += "\""+le.toString()+"\""+",";
					}
					json = json.substring(0, json.length()-1)+"],";
				}
			} else {
				json+="\""+v.toString()+"\",";
			}
		}
		return json.substring(0, json.length()-1)+"}";
	}

	public void fromJSON(JSONObject json) {
		for (Object key: json.keySet()) {
			Object value = json.get(key);
			if (value instanceof JSONArray) {
				List<String> array = new ArrayList<>();
				for (Object obj: (JSONArray)value) {
					array.add((String)obj);
				}
				put((String)key, array);
			} else {
				put((String)key, value);
			}
		}
	}

	private void safePut(String key, Object value) {
		if (value != null)
			put(key, value);
	}
}
