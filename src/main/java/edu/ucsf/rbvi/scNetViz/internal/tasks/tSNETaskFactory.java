package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class tSNETaskFactory extends AbstractTaskFactory {
	final ScNVManager manager;

	public tSNETaskFactory(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public TaskIterator createTaskIterator() {
		return new TaskIterator(new tSNETask(manager));
	}

	public boolean isReady() {
		if (manager.getExperiments() == null || manager.getExperiments().size() == 0)
			return false;
		return true;
	}

}

