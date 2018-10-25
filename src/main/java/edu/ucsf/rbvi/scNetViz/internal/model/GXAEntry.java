package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.ArrayList;
import java.util.List;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

public class GXAEntry {
	String expType;
	String expAccession;
	String expDescription;
	String expDate;
	long expAssays;
	long expContrasts;
	String expSpecies;
	String expKingdom;
	List<String> expFactors;

	public GXAEntry(JSONObject json) {
		expType = (String) json.get("experimentType");
		expAccession = (String) json.get("experimentAccession");
		expDescription = (String) json.get("experimentDescription");
		expDate = (String) json.get("lastUpdate");
		expAssays = (Long) json.get("numberOfAssays");
		expContrasts = (Long) json.get("numberOfContrasts");
		expSpecies = (String) json.get("species");
		expKingdom = (String) json.get("kingdom");
		JSONArray factors = (JSONArray) json.get("experimentalFactors");
		expFactors = new ArrayList<String>();
		for (Object obj: factors) {
			expFactors.add((String)obj);
		}
	}

	public String getType() { return expType; }
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
