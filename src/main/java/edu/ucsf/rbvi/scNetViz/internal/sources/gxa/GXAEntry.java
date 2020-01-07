package edu.ucsf.rbvi.scNetViz.internal.sources.gxa;

import java.util.ArrayList;
import java.util.List;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import edu.ucsf.rbvi.scNetViz.internal.utils.JSONUtils;

public class GXAEntry {
	String expType;
	String expTypeRaw;
	String expAccession;
	String expDescription;
	String expDate;
	long expAssays;
	long expContrasts;
	String expSpecies;
	String expKingdom;
	List<String> expFactors;

	public GXAEntry(JSONObject json) {
		// System.out.println("Entry: "+json.toString());
		expTypeRaw = (String) json.get("rawExperimentType");
		expType = (String) json.get("experimentType");
		expAccession = (String) json.get("experimentAccession");
		expDescription = (String) json.get("experimentDescription");
		expDate = (String) json.get("lastUpdate");
		expAssays = (Long) json.get("numberOfAssays");
		expContrasts = (Long) json.get("numberOfContrasts");
		expSpecies = (String) json.get("species");
		expKingdom = (String) json.get("kingdom");
		JSONArray factors = (JSONArray) json.get("experimentalFactors");
		expFactors = JSONUtils.jsonArrayToList(factors, String.class);
	}

	public String getType() { return expTypeRaw; }
	public String getAccession() { return expAccession; }
	public String getDescription() { return expDescription; }
	public String getDate() { return expDate; }
	public int getAssays() { return (int)expAssays; }
	public int getContrasts() { return (int)expContrasts; }
	public String getSpecies() { return expSpecies; }
	public String getKingdom() { return expKingdom; }
	public List<String> getFactors() { return expFactors; }
	public String toString() {
		return "<html><p style='width: 500px'><b>"+expAccession+"</b>: "+expDescription+"</p></html>";
	}
}
