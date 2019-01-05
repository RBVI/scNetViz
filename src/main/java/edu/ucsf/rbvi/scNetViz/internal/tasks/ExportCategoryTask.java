package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.Arrays;
import java.util.List;
import javax.swing.SwingUtilities;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ContainsTunables;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.json.JSONResult;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class ExportCategoryTask extends AbstractTask {
	final ScNVManager manager;

	@Tunable (description="Experiment accession")
	public ListSingleSelection<Experiment> accession = null;

	public ListSingleSelection<Category> category = null;
	@Tunable (description="Category name",listenForChange="accession")
	public ListSingleSelection<Category> getCategory() {
		if (accession == null) return null;
		List<Category> categories = accession.getSelectedValue().getCategories();
		if (categories != null)
			category = new ListSingleSelection<>(categories);
		else
			category = null;
		return category;
	}
	public void setCategory(ListSingleSelection<Category> category) {
	}

	@ContainsTunables
	public ExportCSVTask exportCSVTask;

	public ExportCategoryTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		accession = new ListSingleSelection<>(manager.getExperiments());
		exportCSVTask = new ExportCSVTask(manager);
	}

	public void run(TaskMonitor monitor) {
		Category cat = category.getSelectedValue();
		Matrix matrix = cat.getMatrix();
		exportCSVTask.setMatrix(matrix);
		insertTasksAfterCurrentTask(exportCSVTask);
	}

	@ProvidesTitle
	public String title() { return "Export Category"; }

}

