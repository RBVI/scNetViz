package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import org.cytoscape.command.AvailableCommands;
import org.cytoscape.command.CommandExecutorTaskFactory;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;

public class ScNVManager {

	final AvailableCommands availableCommands;
	final CommandExecutorTaskFactory ceTaskFactory;
	final TaskManager taskManager;

	final Map<String, Experiment> experimentMap;
	final Map<String, Source> sourceMap;
	final CyServiceRegistrar registrar; 

	public ScNVManager(final CyServiceRegistrar registrar) {
		experimentMap = new HashMap<>();
		sourceMap = new HashMap<>();
		this.registrar = registrar;
		this.availableCommands = registrar.getService(AvailableCommands.class);
		this.ceTaskFactory = registrar.getService(CommandExecutorTaskFactory.class);
		this.taskManager = registrar.getService(TaskManager.class);
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

	public void executeCommand(String namespace, String command, Map<String, Object> args) {
		taskManager.execute(ceTaskFactory.createTaskIterator(namespace, command, args, null));
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

}
