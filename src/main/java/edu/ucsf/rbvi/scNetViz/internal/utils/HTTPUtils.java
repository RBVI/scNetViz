package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.ContentType;
import org.apache.http.entity.StringEntity;
import org.apache.http.entity.mime.MultipartEntityBuilder;
import org.apache.http.entity.mime.content.InputStreamBody;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.apache.http.util.EntityUtils;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class HTTPUtils {
	public static final String WS_URL = "http://webservices.rbvi.ucsf.edu/scnetviz/api/v1/";
	public static final String WS_URL_V2 = "http://webservices.rbvi.ucsf.edu/scnetviz/api/v2/";
	// public static final String WS_URL = "http://localhost:8000/scnetviz/api/v1/";

	public static JSONObject getJSON(String uri, CloseableHttpClient httpclient, 
	                                 TaskMonitor monitor) throws Exception {
		HttpGet httpGet = new HttpGet(uri);
		httpGet.setHeader("Accept", "application/json");
		CloseableHttpResponse response = httpclient.execute(httpGet);
		int statusCode = response.getStatusLine().getStatusCode();
		if (statusCode != 200 && statusCode != 202) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Got "+
			                    response.getStatusLine().getStatusCode()+" code from server");
			return null;
		}

		HttpEntity entity1 = response.getEntity();

		JSONObject jsonResponse = (JSONObject) new JSONParser().parse(new InputStreamReader(entity1.getContent()));
		return jsonResponse;
  }

	public static JSONObject postJSON(String uri, CloseableHttpClient httpclient, 
	                                  String jsonQuery, TaskMonitor monitor) throws Exception {
		HttpPost httpPost = new HttpPost(uri);
		StringEntity entity = new StringEntity(jsonQuery);
		httpPost.setEntity(entity);
		httpPost.setHeader("Accept", "application/json");
		httpPost.setHeader("Content-type", "application/json");

		CloseableHttpResponse response = httpclient.execute(httpPost);
		int statusCode = response.getStatusLine().getStatusCode();
		if (statusCode != 200 && statusCode != 202) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Got "+
			                    response.getStatusLine().getStatusCode()+" code from server");
			return null;
		}

		HttpEntity entity1 = response.getEntity();

		JSONObject jsonResponse = (JSONObject) new JSONParser().parse(new InputStreamReader(entity1.getContent()));
		return jsonResponse;
	}

	public static List<String> postFile(String uri, File f, TaskMonitor monitor) throws Exception {
		System.out.println("postFile: "+uri);
		CloseableHttpClient httpclient = HttpClients.createDefault();
		HttpPost httpPost = new HttpPost(uri);
		MultipartEntityBuilder builder = MultipartEntityBuilder.create();
		builder.addBinaryBody("file", f, ContentType.DEFAULT_BINARY, f.getName());
		HttpEntity entity = builder.build();
		httpPost.setEntity(entity);
		CloseableHttpResponse response = httpclient.execute(httpPost);
		int statusCode = response.getStatusLine().getStatusCode();
		if (statusCode != 200 && statusCode != 202) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Got "+
			                    response.getStatusLine().getStatusCode()+" code from server");
			return null;
		}
		BufferedReader reader = new BufferedReader(new InputStreamReader(response.getEntity().getContent()));
		List<String> strings = new ArrayList<>();
		String line;
		while((line = reader.readLine()) != null) {
			strings.add(line);
		}
		return strings;
	}

	public static List<String> fetchResult(String uri, TaskMonitor monitor) throws Exception {
		CloseableHttpClient httpclient = HttpClients.createDefault();
		HttpGet httpGet = new HttpGet(uri);
		CloseableHttpResponse response1 = httpclient.execute(httpGet);
		int statusCode = response1.getStatusLine().getStatusCode();
		if (statusCode != 200 && statusCode != 202) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Got "+
			                    response1.getStatusLine().getStatusCode()+" code from server");
			return null;
		}
		BufferedReader reader = new BufferedReader(new InputStreamReader(response1.getEntity().getContent()));
		List<String> strings = new ArrayList<>();
		String line;
		while((line = reader.readLine()) != null) {
			strings.add(line);
		}
		return strings;
  }

	public static JSONObject fetchJSON(String uri, TaskMonitor monitor) throws Exception {
		CloseableHttpClient httpclient = HttpClients.createDefault();
		JSONObject obj = fetchJSON(uri, httpclient, monitor);
		httpclient.close();
		return obj;
	}

	public static JSONObject fetchJSON(String uri, CloseableHttpClient httpclient, 
	                                   TaskMonitor monitor) throws Exception {
		System.out.println("Fetching JSON from: "+uri);

		HttpGet httpGet = new HttpGet(uri);
		CloseableHttpResponse response1 = httpclient.execute(httpGet);
		int statusCode = response1.getStatusLine().getStatusCode();
		if (statusCode != 200 && statusCode != 202) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Got "+
			                    response1.getStatusLine().getStatusCode()+" code from server");
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
		return json;
	}

	public static ZipInputStream getZipStream(String uri, 
	                                          CloseableHttpClient httpclient, TaskMonitor monitor) throws Exception {
		System.out.println("Fetching Zip from: "+uri);
		HttpGet httpGet = new HttpGet(uri);
		CloseableHttpResponse response1 = httpclient.execute(httpGet);
		if (response1.getStatusLine().getStatusCode() != 200) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Got "+
			                    response1.getStatusLine().getStatusCode()+" code from server");
			return null;
		}
		ZipInputStream stream = null;
		try {
			HttpEntity entity1 = response1.getEntity();
			stream = new ZipInputStream(entity1.getContent());
		} catch (Exception e) {
			throw e;
		}
		return stream;
	}

	public static String getWebServicesURL(String command, Experiment exp, String args) {
    /*
		String url = WS_URL+command+"?source="+exp.getSource().getName()+"&accession="+
		             exp.getMetadata().get(Metadata.ACCESSION).toString();
    */
		String url = WS_URL_V2+command+"/"+exp.getSource().getName()+"/"+
		             exp.getMetadata().get(Metadata.ACCESSION).toString();
		if (args != null) {
			url += "?"+args;
		}
		return url;
	}

	public static void showWebPage(ScNVManager manager, String uri, TaskMonitor monitor) {
		Map<String, Object> args = new HashMap<>();
		args.put("newTab", "true");
		args.put("id", "GXA");
		args.put("url", uri);

		manager.executeCommand("cybrowser", "dialog", args);

	}
}
