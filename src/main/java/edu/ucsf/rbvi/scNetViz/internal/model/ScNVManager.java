package edu.ucsf.rbvi.scNetViz.internal.model;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.events.SetCurrentNetworkListener;
import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.command.AvailableCommands;
import org.cytoscape.command.CommandExecutorTaskFactory;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.property.CyProperty;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.session.events.SessionAboutToBeSavedEvent;
import org.cytoscape.session.events.SessionAboutToBeSavedListener;
import org.cytoscape.session.events.SessionLoadedEvent;
import org.cytoscape.session.events.SessionLoadedListener;
import org.cytoscape.work.SynchronousTaskManager;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ShowResultsPanelTask;
import edu.ucsf.rbvi.scNetViz.internal.utils.LogUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;
import edu.ucsf.rbvi.scNetViz.internal.view.ScNVCytoPanel;

public class ScNVManager implements SessionAboutToBeSavedListener, SessionLoadedListener {

	final AvailableCommands availableCommands;
	final CommandExecutorTaskFactory ceTaskFactory;
	final TaskManager taskManager;
	final SynchronousTaskManager syncTaskManager;

	final Map<String, Experiment> experimentMap;
	final Map<String, Source> sourceMap;
	final Map<Experiment, ExperimentFrame> frameMap;
	final CyServiceRegistrar registrar; 
	final ScNVSettings settings;
	private ScNVCytoPanel cytoPanel;
	private Icon scNetVizIcon;

	public final static String APP_NAME = "edu.ucsf.rbvi.scNetViz";
	public final static String CATEGORIES = "categories";
	public final static String DIFFEXP = "differential expression";
	public final static String EXPERIMENTS = "Experiments";
	public final static String EXPERIMENTS_FILE = "Experiments.json";
	public final static String SOURCE_NAME = "source name";

	public ScNVManager(final CyServiceRegistrar registrar) {
		experimentMap = new HashMap<>();
		sourceMap = new HashMap<>();
		frameMap = new HashMap<>();
		this.registrar = registrar;
		this.availableCommands = registrar.getService(AvailableCommands.class);
		this.ceTaskFactory = registrar.getService(CommandExecutorTaskFactory.class);
		this.taskManager = registrar.getService(TaskManager.class);
		this.syncTaskManager = registrar.getService(SynchronousTaskManager.class);
		settings = new ScNVSettings();
		scNetVizIcon = new ImageIcon(getClass().getResource("/images/scNetViz.png"));

		registrar.registerService(this, SessionAboutToBeSavedListener.class, new Properties());
		registrar.registerService(this, SessionLoadedListener.class, new Properties());
	}

	public void addSource(Source source) {
		sourceMap.put(source.getName(), source);
	}

	public void addSource(String name, Source source) {
		sourceMap.put(name, source);
	}

	public List<Source> getSources() {
		return new ArrayList<>(sourceMap.values());
	}

	public Source getSource(String name) {
		return sourceMap.get(name);
	}

	public void addExperiment(String accession, Experiment exp) {
		experimentMap.put(accession, exp);
	}

	public void deleteExperiment(String accession) {
		if (experimentMap.containsKey(accession)) {
			Experiment exp = experimentMap.get(accession);
			experimentMap.remove(accession);
			if (frameMap.containsKey(exp)) {
				frameMap.get(exp).dispose();
				frameMap.remove(exp);
			}

			//TODO: what about results panel?
			if (cytoPanel != null) {
				if (cytoPanel.getExperiment().equals(exp)) {
					unregisterService(cytoPanel, CytoPanelComponent.class);
					cytoPanel = null;
				}
			}
		}
	}

	public Experiment getExperiment(String accession) {
		if (experimentMap.containsKey(accession)) return experimentMap.get(accession);
		return null;
	}

	public List<Experiment> getExperiments() {
		return new ArrayList<>(experimentMap.values());
	}

	public Set<String> getExperimentAccessions() {
		return experimentMap.keySet();
	}

	public void addExperimentFrame(Experiment experiment, ExperimentFrame expFrame) {
		frameMap.put(experiment, expFrame);
	}

	public ExperimentFrame getExperimentFrame(Experiment experiment) {
		return frameMap.get(experiment);
	}

	public String getSetting(ScNVSettings.SETTING setting) {
		return settings.getSetting(setting);
	}

	public void setSetting(ScNVSettings.SETTING setting, double value) {
		setSetting(setting, String.valueOf(value));
	}

	public void setSetting(ScNVSettings.SETTING setting, int value) {
		setSetting(setting, String.valueOf(value));
	}

	public void setSetting(ScNVSettings.SETTING setting, boolean value) {
		setSetting(setting, String.valueOf(value));
	}

	public void setSetting(ScNVSettings.SETTING setting, String value) {
		settings.setSetting(setting, value);
	}

	public void setCytoPanel(ScNVCytoPanel panel) {
		this.cytoPanel = panel;
	}

	public Icon getIcon() { return scNetVizIcon; }

	public ScNVCytoPanel getCytoPanel() { return this.cytoPanel; }

	public void executeCommand(String namespace, String command, Map<String, Object> args, boolean synchronous) {
		executeCommand(namespace, command, args, null, synchronous);
	}

	public void executeCommand(String namespace, String command, Map<String, Object> args) {
		executeCommand(namespace, command, args, null, false);
	}

	public void executeCommand(String namespace, String command, Map<String, Object> args, 
	                           TaskObserver observer, boolean synchronous) {
		List<String> commands = availableCommands.getCommands(namespace);
		if (!commands.contains(command)) {
			LogUtils.warn("Command "+namespace+" "+command+" isn't available");
			return;
		}

		if (synchronous)
			syncTaskManager.execute(ceTaskFactory.createTaskIterator(namespace, command, args, observer));
		else
			taskManager.execute(ceTaskFactory.createTaskIterator(namespace, command, args, observer));
	}

	public boolean haveNamespace(String namespace) {
		List<String> namespaces = availableCommands.getNamespaces();
		if (namespaces.contains(namespace)) return true;
		return false;
	}

	public boolean haveCommand(String namespace, String command) {
		if (!haveNamespace(namespace)) return false;
		List<String> commands = availableCommands.getCommands(namespace);
		if (commands.contains(command)) return true;
		return false;
	}

	public void executeTasks(TaskIterator tasks) {
		taskManager.execute(tasks);
	}

	public void executeTasks(TaskIterator tasks, TaskObserver observer) {
		taskManager.execute(tasks, observer);
	}

	public void executeTasks(TaskFactory factory) {
		taskManager.execute(factory.createTaskIterator());
	}

	public void executeTasks(TaskFactory factory, TaskObserver observer) {
		taskManager.execute(factory.createTaskIterator(), observer);
	}

	public <S> S getService(Class<S> serviceClass) {
		return registrar.getService(serviceClass);
	}

	public <S> S getService(Class<S> serviceClass, String filter) {
		return registrar.getService(serviceClass, filter);
	}

	public void registerService(Object service, Class<?> serviceClass, Properties props) {
		registrar.registerService(service, serviceClass, props);
	}

	public void unregisterService(Object service, Class<?> serviceClass) {
		registrar.unregisterService(service, serviceClass);
	}

	// See if we have data in the session, and load it if we do
	public void handleEvent(SessionLoadedEvent e) {
		// First, if we have a results panel, unregister it
		if (cytoPanel != null) {
			unregisterService(cytoPanel, CytoPanelComponent.class);
			unregisterService(cytoPanel, SetCurrentNetworkListener.class);
			cytoPanel = null;
		}
		Map<String,List<File>> appFiles = e.getLoadedSession().getAppFileListMap();
		if (!appFiles.containsKey(APP_NAME))
			return;

		List<File> scNvFiles = appFiles.get(APP_NAME);
		Map<String, File> fileMap = new HashMap<>();
		for (File f: scNvFiles) {
			System.out.println("File map has file: "+f.getName());
			fileMap.put(f.getName(),f);
		}

		if (!fileMap.containsKey(EXPERIMENTS_FILE))
			return;

		JSONParser parser = new JSONParser();
		JSONObject jsonExperiment;
		try {
			jsonExperiment = (JSONObject) parser.parse(new FileReader(fileMap.get(EXPERIMENTS_FILE)));
		} catch(Exception ioe) {
			return;
		}

		JSONArray experiments = (JSONArray) jsonExperiment.get(EXPERIMENTS);
		for (Object exp: experiments) {
			JSONObject jsonExp = (JSONObject) exp;
			Experiment experiment = getExperimentFromSession(jsonExp, fileMap);
			// System.out.println("Loaded expermient: "+experiment);
			if (experiment == null) continue;

			if (jsonExp.containsKey(CATEGORIES)) {
				// System.out.println("Getting categories");
				JSONArray categories = (JSONArray) jsonExp.get(CATEGORIES);
				for (Object cat: categories) {
					Category category = getCategoryFromSession((JSONObject)cat, experiment, fileMap);
				}
			}

			if (jsonExp.containsKey(DIFFEXP)) {
				getDiffExpFromSession((JSONObject)jsonExp.get(DIFFEXP), experiment, fileMap);
			}
		}

		CyNetwork network = getService(CyApplicationManager.class).getCurrentNetwork();
		if (network != null) {
			Experiment exp = ModelUtils.getExperimentFromNetwork(this, network);
			if (exp != null) {
				// Now, show the results panel
				Task resultsPanelTask = new ShowResultsPanelTask(this, exp);
				executeTasks(new TaskIterator(resultsPanelTask));
			}
		}
	}

	// We need to save all of our experiment, category, and DE tables
	public void handleEvent(SessionAboutToBeSavedEvent e) {
		String tmpDir = System.getProperty("java.io.tmpdir");
		File jsonFile = new File(tmpDir, EXPERIMENTS_FILE);

		try {
			FileOutputStream fos = new FileOutputStream(jsonFile);
			OutputStreamWriter osw = new OutputStreamWriter(fos, "utf-8");
			BufferedWriter writer = new BufferedWriter(osw);

			writer.write("{\""+EXPERIMENTS+"\":[");
			List<File> files = new ArrayList<File>();
			int expNumber = experimentMap.keySet().size();
			for (String accession: experimentMap.keySet()) {
				Experiment experiment = experimentMap.get(accession);
				writer.write(experiment.toJSON());
				if (expNumber-- > 1)
					writer.write(",\n");
				try {
					experiment.createSessionFiles(accession, files);
				} catch (Exception create) {
				}
			}
			writer.write("]}\n");
			writer.close();
			osw.close();
			fos.close();
			files.add(jsonFile);

			try {
				e.addAppFiles(APP_NAME, files);
			} catch (Exception add) {
			}
		} catch (Exception jsonException) {
		}
	}

	private Experiment getExperimentFromSession(JSONObject jsonExp, Map<String,File> fileMap) {
		String src = (String)jsonExp.get(SOURCE_NAME);
		System.out.println("Getting an experiment from "+src);
		if (sourceMap.containsKey(src))
			return sourceMap.get(src).loadExperimentFromSession(jsonExp, fileMap);
		return null;
	}

	private Category getCategoryFromSession(JSONObject jsonCategory, Experiment experiment, Map<String,File> fileMap) {
		String src = (String)jsonCategory.get(SOURCE_NAME);
		System.out.println("Getting a category from "+src);
		if (sourceMap.containsKey(src))
			return sourceMap.get(src).loadCategoryFromSession(jsonCategory, experiment, fileMap);
		return null;
	}

	private void getDiffExpFromSession(JSONObject jsonDiffExp, Experiment experiment, Map<String,File> fileMap) {
		System.out.println("Getting differential expression for experiment "+experiment);
		experiment.getSource().loadDiffExpFromSession(jsonDiffExp, experiment, fileMap);
	}
}
