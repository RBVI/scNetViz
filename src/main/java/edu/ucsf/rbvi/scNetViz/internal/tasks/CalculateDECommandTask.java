package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings;
import edu.ucsf.rbvi.scNetViz.internal.view.DiffExpTab;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;

// Tunable to choose experiment?

public class CalculateDECommandTask extends AbstractTask implements ObservableTask {
	final ScNVManager manager;
	DifferentialExpression diffExp;

	// Tunables:
	//    Experiment accession
	//    Category category
	//    int categoryRow
	//    double minPctCutoff
	//    double logGERCutoff
	@Tunable (description="Experiment accession")
	public ListSingleSelection<Experiment> accession = null;

	public ListSingleSelection<Category> category = null;
	@Tunable (description="Category name",listenForChange="accession")
	public ListSingleSelection<Category> getCategory() {
		if (accession == null) return null;
		List<Category> categories = accession.getSelectedValue().getCategories();
		if (categories != null) {
			category = new ListSingleSelection<>(categories);
			category.setSelectedValue(accession.getSelectedValue().getDefaultCategory());
		} else
			category = null;
		return category;
	}
	public void setCategory(ListSingleSelection<Category> category) {
	}

	public int categoryRow = -1;
	@Tunable (description="The category row to use", listenForChange="category")
	public int getCategoryRow() {
		if (categoryRow == -1 && category != null) {
			Category cat = category.getSelectedValue();
			categoryRow = cat.getSelectedRow();
			if (categoryRow < 0)
				categoryRow = cat.getDefaultRow();
			if (categoryRow >= 0)
				cat.setSelectedRow(categoryRow);
		}
		return categoryRow;
	}
	public void setCategoryRow(int value) {
		categoryRow = value;
	}

	@Tunable (description="Min.pct cutoff")
	public double minPctCutoff;

	@Tunable (description="logGER cutoff")
	public double logGERCutoff;

	public CalculateDECommandTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		accession = new ListSingleSelection<>(manager.getExperiments());
		logGERCutoff = Double.parseDouble(manager.getSetting(ScNVSettings.SETTING.DE_FC_CUTOFF));
		minPctCutoff = Double.parseDouble(manager.getSetting(ScNVSettings.SETTING.DE_MIN_PCT_CUTOFF));
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("Calculating Differential Expression");
		Category cat;
		Experiment exp;

		if (manager.getExperiments().size() == 0) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "No experiments loaded");
			throw new RuntimeException("No experiments loaded");
		}

		exp = accession.getSelectedValue();

		if (category == null) {
			cat = exp.getDefaultCategory();
			if (cat == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "No category specified and no default is defined");
				throw new RuntimeException("No category specified and no default is defined");
			}
		} else
			cat = category.getSelectedValue();
			
		if (categoryRow == -1) {
			categoryRow = cat.getSelectedRow();
			if (categoryRow < 0)
				categoryRow = cat.getDefaultRow();
		}

    if (categoryRow == -1) {
		  monitor.showMessage(TaskMonitor.Level.ERROR, 
		                      "Unable to complete differential expression calculation: no default row available.");
      return;
    }

		try {
			/*
			System.out.println("experiment = "+accession.getSelectedValue());
			System.out.println("category = "+category.getSelectedValue());
			System.out.println("categoryRow = "+categoryRow);
			System.out.println("minPctCutoff = "+minPctCutoff);
			System.out.println("logGERCutoff = "+logGERCutoff);
			*/
			diffExp = new DifferentialExpression(manager, cat, 
			                                     categoryRow, minPctCutoff/100, logGERCutoff);
			// System.out.println("diffExp = "+diffExp);
			exp.setDiffExp(diffExp);
			ExperimentFrame frame = manager.getExperimentFrame(exp);
			DiffExpTab diffETab = new DiffExpTab(manager, exp, frame, cat, diffExp);
			frame.addDiffExpContent("Diff Exp", diffETab);
			monitor.showMessage(TaskMonitor.Level.INFO, 
			                    "Calculations complete");
		} catch (Exception ex) {
			ex.printStackTrace();
			monitor.showMessage(TaskMonitor.Level.ERROR, 
			                    "Unable to complete differential expression calculation: "+ex.getMessage());
		}
	}

	public <R> R getResults(Class<? extends R> clazz) {
		if (clazz.equals(DifferentialExpression.class))
			return (R)diffExp;
		return null;
	}
}
