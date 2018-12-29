package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class SelectTaskFactory extends AbstractTaskFactory {
	final ScNVManager manager;

	public SelectTaskFactory(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public TaskIterator createTaskIterator() {
		return new TaskIterator(new SelectTask(manager));
	}

	public boolean isReady() {
		if (manager.getExperiments() == null || manager.getExperiments().size() == 0)
			return false;
		return true;
	}

}

