package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

// Tunable to choose experiment?

public class CalculateDETask extends AbstractTask implements ObservableTask {
	final ScNVManager manager;

	// FIXME: these should be Tunables at some point
	Category category = null;
	double dDRCutoff;
	double log2FCCutoff;
	DifferentialExpression diffExp = null;

	public CalculateDETask(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public CalculateDETask(final ScNVManager manager, Category category, double dDRCutoff, double log2FCCutoff) {
		super();
		this.manager = manager;
		this.category = category;
		this.dDRCutoff = dDRCutoff;
		this.log2FCCutoff = log2FCCutoff;
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("Calculating Differential Expression");
		int row = category.getSelectedRow();
		System.out.println("Calculating differential expression: row = "+row);
		if (row < 0) {
			row = category.getDefaultRow();
			if (row < 0) {
				// Nothing selected and no default
				monitor.showMessage(TaskMonitor.Level.ERROR, 
			                      "No row was selected and this category does not have a default row");
				return;
			}
			category.setSelectedRow(row);
		}
		try {
			diffExp = new DifferentialExpression(manager, category, row, dDRCutoff, log2FCCutoff);
			category.getExperiment().setDiffExp(diffExp);
		} catch (Exception exp) {
			exp.printStackTrace();
			monitor.showMessage(TaskMonitor.Level.ERROR, 
			                    "Unable complete differential expression calculation: "+exp.getMessage());
		}
	}

	public <R> R getResults(Class<? extends R> clazz) {
		if (clazz.equals(DifferentialExpression.class))
			return (R)diffExp;
		return null;
	}
}
