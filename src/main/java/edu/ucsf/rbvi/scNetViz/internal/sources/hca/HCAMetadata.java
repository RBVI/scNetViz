package edu.ucsf.rbvi.scNetViz.internal.sources.hca;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
	public static String FULL_DESCRIPTION = "name";
	public static String MATRIX = "matrix";
	public static String SELECTED_CELL_TYPES = "cellTypes";

	public HCAMetadata(String organ, JSONObject json) {
		super();
		put(Metadata.ACCESSION, (String)json.get("entryId")+"-"+organ);
		getProjects((JSONArray)json.get("projects"), organ);
		getSuspensions((JSONArray)json.get("cellSuspensions"), organ);
  }

	private void getProjects(JSONArray projects, String organ) {
		JSONObject project = (JSONObject)projects.get(0); // Only dealing with one project at a time, here
		put(DESCRIPTION, (String)project.get("projectTitle"));
		put(SHORT_NAME, (String)project.get("projectShortname"));
		put(FULL_DESCRIPTION, (String)project.get("projectDescription"));
		put(GEO, getArray((JSONArray)project.get("geoSeriesAccessions")));
		put(ARRAY_EXPRESS, getArray((JSONArray)project.get("arrayExpressAccessions")));
		put(INSDC, getArray((JSONArray)project.get("insdcProjectAccessions")));
		put(LABORATORY, getArray((JSONArray)project.get("laboratory")));
    getMatrix((JSONObject)project.get("matrices"), organ);
  }

	private void getSuspensions(JSONArray suspensions, String organ) {
    for (Object susp: suspensions) {
      JSONObject suspension = (JSONObject) susp;
      JSONArray organs = (JSONArray)suspension.get("organ");
      for (Object org: organs) {
        if (organ.equals((String)org)) {
		      put(ASSAYS, ((Long)suspension.get("totalCells")).longValue());
          put(SELECTED_CELL_TYPES, getArray((JSONArray)suspension.get("selectedCellType")));
          put(ORGANS, Collections.singletonList(organ));
        }
      }
    }
    return;
  }


	private static List<String> getArray(JSONArray jsonArray) {
		List<String> array = new ArrayList<String>();
		for (Object obj: jsonArray) {
			array.add((String)obj);
		}
		return array;
	}

  private void getMatrix(JSONObject jsonMatrices, String organ) {
    JSONObject obj = (JSONObject) jsonMatrices.get("genusSpecies");
    for (Object spObj: obj.keySet()) {
      String species = (String)spObj;
      JSONObject xobj = (JSONObject)obj.get(spObj);
      if (xobj.containsKey("organ")) {
        JSONObject organs = (JSONObject)xobj.get("organ");
        if (organs.containsKey(organ)) {
          JSONObject orgObj = (JSONObject)organs.get(organ);
          JSONObject lcaObj = (JSONObject)orgObj.get("libraryConstructionApproach");
          for (Object lcaType: lcaObj.keySet()) {
            JSONArray matrices = (JSONArray)lcaObj.get(lcaType);
            for (Object matObj: matrices) {
              String name = (String)((JSONObject)matObj).get("name");
              if (name.contains("mtx.zip")) {
                put(Metadata.TYPE, (String)lcaType);
                put(MATRIX, (String)((JSONObject)matObj).get("url"));
                put(Metadata.SPECIES, species);
                return;
              }
            }
          }
        }
      } else if (xobj.containsKey("libraryConstructionApproach")) {
        JSONObject lcas = (JSONObject)xobj.get("libraryConstructionApproach");
        for (Object lca: lcas.keySet()) {
          JSONObject lcaObj = (JSONObject)lcas.get(lca);
          if (!lcaObj.containsKey("organ"))
            continue;
          JSONObject organs = (JSONObject) lcaObj.get("organ");
          if (!organs.containsKey(organ))
            continue;
          JSONArray matrices = (JSONArray)organs.get(organ);
          for (Object matObj: matrices) {
            String name = (String)((JSONObject)matObj).get("name");
            if (name.contains("mtx.zip")) {
              put(Metadata.TYPE, (String)lca);
              put(Metadata.SPECIES, species);
              put(MATRIX, (String)((JSONObject)matObj).get("url"));
              return;
            }
          }
        }
      }
    }
    return;
  }

  public static List<String> getOrgans(JSONObject hit) {
		JSONArray suspensions = (JSONArray)hit.get("cellSuspensions");
    List<String> organList = new ArrayList<>();
    for (Object suspension: suspensions) {
      JSONObject susp = (JSONObject) suspension;
		  organList.addAll(getArray((JSONArray)susp.get("organ")));
    }
    return organList;
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
