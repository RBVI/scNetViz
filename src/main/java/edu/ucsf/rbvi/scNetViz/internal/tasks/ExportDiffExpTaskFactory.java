package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class ExportDiffExpTaskFactory extends AbstractTaskFactory {
	final ScNVManager manager;

	public ExportDiffExpTaskFactory(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public TaskIterator createTaskIterator() {
		return new TaskIterator(new ExportDiffExpTask(manager));
	}

	public boolean isReady() {
		return true;
	}

}

