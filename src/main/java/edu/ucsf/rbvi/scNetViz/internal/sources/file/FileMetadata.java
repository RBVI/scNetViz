package edu.ucsf.rbvi.scNetViz.internal.sources.file;

import java.io.File;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;

public class FileMetadata extends HashMap<String, Object> implements Metadata {
	public static String ASSAYS = "assays";
	public static String FACTORS = "factors";
	public static String FILE = "file";

	public FileMetadata(File file) {
		super();
		put(ACCESSION, file.getName());
		put(TYPE,"File");
		put(DESCRIPTION, "Experiment loaded from file");
		put(FILE, file);
	}

  public FileMetadata() {
    super();
  }

	public String toString() {
		return get(ACCESSION)+": "+get(DESCRIPTION);
	}

	public String toHTML() {
		return "<html><p style='width: 500px'><b>"+get(ACCESSION)+"</b>: "+get(DESCRIPTION)+"</p></html>";
	}

	public String toJSON() {
		String json = "{";
		for (String key: keySet()) {
			Object v = get(key);
			json += "\""+key+"\":";
			if (v instanceof List) {
				json += "[";
				for (Object lv: (List)v) {
					json += "\""+lv.toString()+"\""+",";
				}
				json = json.substring(0, json.length()-1)+"],";
			} else {
				json+="\""+v.toString()+"\",";
			}
		}
		return json.substring(0, json.length()-1)+"}";
	}

  /*
  Getting experiment: {"file":"\/var\/tmp\/TCGA_BRCA\/TCGA_BRCA.expression.fpkm.zip","species":"Homo sapiens","description":"Experiment loaded from file","accession":"TCGA_BRCA.expression.fpkm.zip","type":"File"}
  */
	public static FileMetadata fromJSON(JSONObject metadata) {
    File f = new File((String)metadata.get("file"));
    FileMetadata fm = new FileMetadata(f);
    return fm;
  }
}
