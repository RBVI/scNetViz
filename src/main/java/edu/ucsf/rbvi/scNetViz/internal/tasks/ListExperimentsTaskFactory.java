package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class ListExperimentsTaskFactory extends AbstractTaskFactory {
	final ScNVManager manager;

	public ListExperimentsTaskFactory(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public TaskIterator createTaskIterator() {
		return new TaskIterator(new ListExperimentsTask(manager));
	}

	public boolean isReady() {
		return true;
	}

}

