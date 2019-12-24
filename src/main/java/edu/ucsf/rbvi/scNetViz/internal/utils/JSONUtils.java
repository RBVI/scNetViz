package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.util.ArrayList;
import java.util.List;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

public class JSONUtils {
	public static <T> List<T> jsonArrayToList(JSONArray array, Class<T> type) {
		List<T> returnList = new ArrayList<>();

		for (Object obj: array) {
			returnList.add((T)obj);
		}
		return returnList;
	}
}
