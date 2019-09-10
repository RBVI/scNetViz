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
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;

// Tunable to choose experiment?

public class ProcessAllTask extends AbstractTask implements TaskObserver {
	final ScNVManager manager;
	Experiment experiment = null;
	// int maxGenes = 100000;
	int maxGenes = 50;

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

		double logFC = 0.5;
		if (haveSetting(SETTING.DE_FC_CUTOFF))
			logFC =	Double.parseDouble(manager.getSetting(SETTING.DE_FC_CUTOFF));
		double dDRThreshold = 0.1;
		if (haveSetting(SETTING.DE_MIN_PCT_CUTOFF))
			dDRThreshold = Double.parseDouble(manager.getSetting(SETTING.DE_MIN_PCT_CUTOFF))/100.0;

		TaskIterator ti = new TaskIterator(new CalculateDETask(manager, defaultCategory, dDRThreshold, logFC));
		manager.executeTasks(ti, this);
	}

	private boolean haveSetting(SETTING setting) {
		if (manager.getSetting(setting) != null && manager.getSetting(setting).length() > 0)
			return true;
		return false;
	}

	@Override
	public void allFinished(FinishStatus status) {
	}

	@Override
	public void taskFinished(ObservableTask obsTask) {
		if (obsTask instanceof CalculateDETask) {
			DifferentialExpression diffExp = obsTask.getResults(DifferentialExpression.class);
			double pValue = Double.parseDouble(manager.getSetting(SETTING.NET_PV_CUTOFF));
			double log2FCCutoff = Double.parseDouble(manager.getSetting(SETTING.NET_FC_CUTOFF));
			int maxGenes = Integer.parseInt(manager.getSetting(SETTING.MAX_GENES));
			boolean positiveOnly = Boolean.parseBoolean(manager.getSetting(SETTING.POSITIVE_ONLY));
			// TODO: use a reasonable default for maxGenes: 500?
			TaskIterator ti = new TaskIterator(new CreateNetworkTask(manager, diffExp, pValue, log2FCCutoff, maxGenes, positiveOnly, maxGenes));
			manager.executeTasks(ti);
		}
	}
}
