package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class ExportCSVTaskFactory extends AbstractTaskFactory {
	final ScNVManager manager;
	final Matrix matrix;

	public ExportCSVTaskFactory(final ScNVManager manager, final Matrix matrix) {
		super();
		this.manager = manager;
		this.matrix = matrix;
	}

	public TaskIterator createTaskIterator() {
		return new TaskIterator(new ExportCSVTask(manager, matrix));
	}

	public boolean isReady() {
		if (manager.getExperiments() == null || manager.getExperiments().size() == 0)
			return false;
		return true;
	}

}

