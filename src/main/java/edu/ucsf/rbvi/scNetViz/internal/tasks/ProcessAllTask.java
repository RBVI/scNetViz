package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

// Tunable to choose experiment?

public class ProcessAllTask extends AbstractTask {
	final ScNVManager manager;
	Experiment experiment = null;

	public ProcessAllTask(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public ProcessAllTask(final ScNVManager manager, Experiment experiment) {
		super();
		this.manager = manager;
		this.experiment = experiment;
	}

	public void run(TaskMonitor monitor) {
	}
}
