package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

// Tunable to choose experiment?

public class ProcessAllTask extends AbstractTask implements TaskObserver {
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
		Category defaultCategory = experiment.getDefaultCategory();
		if (defaultCategory == null) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "No default category for this experiment");
			return;
		}

		int defaultRow = defaultCategory.getDefaultRow();
		if (defaultRow == -1) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "No default row for this category");
			return;
		}

		defaultCategory.setSelectedRow(defaultRow);

		TaskIterator ti = new TaskIterator(new CalculateDETask(manager, defaultCategory, 0.1, 0.5));
		manager.executeTasks(ti, this);
	}

	@Override
	public void allFinished(FinishStatus status) {
	}

	@Override
	public void taskFinished(ObservableTask obsTask) {
		if (obsTask instanceof CalculateDETask) {
			DifferentialExpression diffExp = obsTask.getResults(DifferentialExpression.class);
			TaskIterator ti = new TaskIterator(new CreateNetworkTask(manager, diffExp, 0.05, 1.0, -1));
			manager.executeTasks(ti);
		}
	}
}
