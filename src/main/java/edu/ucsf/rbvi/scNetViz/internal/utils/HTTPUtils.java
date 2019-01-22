package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.io.BufferedReader;
import java.io.InputStreamReader;

import java.util.HashMap;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class HTTPUtils {
	public static JSONObject fetchJSON(String uri, TaskMonitor monitor) {
		System.out.println("Fetching JSON from: "+uri);
		try {
			CloseableHttpClient httpclient = HttpClients.createDefault();
			HttpGet httpGet = new HttpGet(uri);
			CloseableHttpResponse response1 = httpclient.execute(httpGet);
			if (response1.getStatusLine().getStatusCode() != 200) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "Got "+response1.getStatusLine().getStatusCode()+" code from server");
				return null;
			}
			HttpEntity entity1 = response1.getEntity();
			BufferedReader reader = new BufferedReader(new InputStreamReader(entity1.getContent()));
			/*
			String line;
			while ((line = reader.readLine()) != null) {
				System.out.println(line);
			}
			*/
			JSONObject json = (JSONObject) new JSONParser().parse(reader);
			httpclient.close();
			return json;
		} catch (Exception e) {}
		return null;
	}

	public static void showWebPage(ScNVManager manager, String uri, TaskMonitor monitor) {
		Map<String, Object> args = new HashMap<>();
		args.put("newTab", "true");
		args.put("id", "GXA");
		args.put("url", uri);

		manager.executeCommand("cybrowser", "dialog", args);

	}
}
