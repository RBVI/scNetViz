package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.Arrays;
import java.util.List;
import javax.swing.SwingUtilities;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.json.JSONResult;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class GetExperimentTask extends AbstractTask implements ObservableTask {
	final ScNVManager manager;
	Experiment experiment = null;

	@Tunable (description="Experiment to show")
	public ListSingleSelection<Experiment> accession = null;

	public GetExperimentTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		accession = new ListSingleSelection<>(manager.getExperiments());
	}

	public void run(TaskMonitor monitor) {
		experiment = accession.getSelectedValue();
	}

	@ProvidesTitle
	public String title() { return "GetExperiments"; }

	@Override
	public List<Class<?>> getResultClasses() {
		return Arrays.asList(String.class, JSONResult.class, Experiment.class);
	}

	@Override
	public <R> R getResults(Class<? extends R> type) {
		if (type.equals(Experiment.class))
			return (R)experiment;
		else if (type.equals(JSONResult.class)) {
			JSONResult res = () -> {
				return experiment.toJSON();
			};
			return (R)res;
		} else {
			StringBuilder builder = new StringBuilder();
			builder.append("source: "+experiment.getSource().toString()+"\n");
			builder.append("accession: "+experiment.getMetadata().get(Metadata.ACCESSION).toString()+"\n");
			builder.append("species: "+experiment.getSpecies().toString()+"\n");
			builder.append("description: "+experiment.getMetadata().get(Metadata.DESCRIPTION).toString()+"\n");
			builder.append("rows: "+experiment.getMatrix().getNRows()+"\n");
			builder.append("columns: "+experiment.getMatrix().getNCols()+"\n");
			List<Category> categories = experiment.getCategories();
			builder.append("categories: ");
			for (Category cat: categories) {
				builder.append("    "+cat.toString()+"\n");
			}
			return (R)builder.toString();
		}
	}
}

