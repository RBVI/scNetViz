package edu.ucsf.rbvi.scNetViz.internal.sources.hca;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;

// Format of the JSON:
//  'entryId': 'cddab57b-6868-4be4-806f-395ed9dd635a', 
//  'samples': [{'sampleEntityType': ['specimens'], 'id': ['DID_scRSq05_pancreas', 'DID_scRSq07_pancreas', 'DID_scRSq04_pancreas', 'DID_scRSq08_pancreas', 'DID_scRSq06_pancreas', 'DID_scRSq01_pancreas', 'DID_scRSq02_pancreas', 'DID_scRSq03_pancreas'], 'organ': ['pancreas'], 'organPart': ['islet of Langerhans'], 'disease': ['normal'], 'preservationMethod': [None], 'source': ['specimen_from_organism']}], 
//  'specimens': [{'id': ['DID_scRSq05_pancreas', 'DID_scRSq07_pancreas', 'DID_scRSq04_pancreas', 'DID_scRSq08_pancreas', 'DID_scRSq06_pancreas', 'DID_scRSq01_pancreas', 'DID_scRSq02_pancreas', 'DID_scRSq03_pancreas'], 'organ': ['pancreas'], 'organPart': ['islet of Langerhans'], 'disease': ['normal'], 'preservationMethod': [None], 'source': ['specimen_from_organism']}], 
//  'cellLines': [], 
//  'donorOrganisms': [{'id': ['DID_scRSq05', 'DID_scRSq07', 'DID_scRSq08', 'DID_scRSq06', 'DID_scRSq02', 'DID_scRSq01', 'DID_scRSq03', 'DID_scRSq04'], 'genusSpecies': ['Homo sapiens'], 'organismAge': ['6', '1', '44', '38', '5', '54', '22', '21'], 'organismAgeUnit': ['month', 'year'], 'biologicalSex': ['male', 'female'], 'disease': None}], 
//  'organoids': [], 
//  'cellSuspensions': [{'organ': ['pancreas'], 'organPart': ['islet of Langerhans'], 'selectedCellType': [], 'totalCells': 2544}], 
//  'fileTypeSummaries': [{'fileType': 'fastq.gz', 'count': 5088, 'totalSize': 359718941993}, {'fileType': 'txt', 'count': 12720, 'totalSize': 47079216}, {'fileType': 'csv', 'count': 17808, 'totalSize': 237750043}, {'fileType': 'bam', 'count': 5088, 'totalSize': 1609285774051}, {'fileType': 'bai', 'count': 2544, 'totalSize': 4971850256}, {'fileType': 'results', 'count': 5088, 'totalSize': 67331466554}, {'fileType': 'matrix', 'count': 2544, 'totalSize': 376512}], 
//           }
//

public class HCAMetadata extends HashMap<String, Object> implements Metadata {
	public static String DONOR_COUNT = "donorCount";
	public static String DISEASE = "disease";
	public static String ORGANS = "organs";
	public static String ASSAYS = "assays";
	public static String GEO = "geoAccession";
	public static String ARRAY_EXPRESS = "arrayExpress";
	public static String INSDC = "inscd";
	public static String LABORATORY = "laboratory";
	public static String SHORT_NAME = "name";

	public HCAMetadata(JSONObject json) {
		super();

		put(Metadata.ACCESSION, (String)json.get("entryId"));
		getProtocols((JSONArray)json.get("protocols"));
		getProjects((JSONArray)json.get("projects"));
		getSummary((JSONObject)json.get("projectSummary"));
	}

	private void getProtocols(JSONArray protocols) {
		// {'protocols': [{'libraryConstructionApproach': ['Smart-seq2'], 'instrumentManufacturerModel': ['Illumina NextSeq 500'], 'pairedEnd': [True], 'workflow': ['smartseq2_v2.3.0', 'smartseq2_v2.4.0'], 'assayType': []}], 
		JSONObject protocol = (JSONObject)protocols.get(0);
	}

	private void getProjects(JSONArray projects) {
		//  'projects': [{'projectTitle': 'Single cell transcriptome analysis of human pancreas reveals transcriptional signatures of aging and somatic mutation patterns.', 'projectShortname': 'Single cell transcriptome analysis of human pancreas', 'laboratory': ['Human Cell Atlas Data Coordination Platform', 'Molecular Atlas'], 'arrayExpressAccessions': [], 'geoSeriesAccessions': ['GSE81547'], 'insdcProjectAccessions': ['SRP075496'], 'insdcStudyAccessions': []}], 
		JSONObject project = (JSONObject)projects.get(0);
		put(DESCRIPTION, (String)project.get("projectTitle"));
		put(SHORT_NAME, (String)project.get("projectShortname"));
		put(GEO, getArray((JSONArray)project.get("geoSeriesAccessions")));
		put(ARRAY_EXPRESS, getArray((JSONArray)project.get("arrayExpressAccessions")));
		put(INSDC, getArray((JSONArray)project.get("insdcProjectAccessions")));
		put(LABORATORY, getArray((JSONArray)project.get("laboratory")));
	}

	private void getSummary(JSONObject summary) {
		//  'projectSummary': {'donorCount': 8, 'totalCellCount': 2544.0, 'organTypes': ['pancreas'], 'cellCountSummaries': [{'organType': ['pancreas'], 'countOfDocsWithOrganType': 1, 'totalCellCountByOrgan': 2544.0}], 'genusSpecies': ['Homo sapiens'], 'libraryConstructionApproach': ['Smart-seq2'], 'disease': ['normal']}
		put(ASSAYS, ((Double)summary.get("totalCellCount")).longValue());
		put(Metadata.SPECIES, (String)((JSONArray)summary.get("genusSpecies")).get(0));
		put(Metadata.TYPE, (String)((JSONArray)summary.get("libraryConstructionApproach")).get(0));
		put(ORGANS, getArray((JSONArray)summary.get("organTypes")));
	}

	private List<String> getArray(JSONArray jsonArray) {
		List<String> array = new ArrayList<String>();
		for (Object obj: jsonArray) {
			array.add((String)obj);
		}
		return array;
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
}
