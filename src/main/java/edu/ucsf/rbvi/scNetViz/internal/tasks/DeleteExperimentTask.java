package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
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

public class DeleteExperimentTask extends AbstractTask implements ObservableTask {
	final ScNVManager manager;
	String experiment = null;

	@Tunable (description="Experiment to delete")
	public ListSingleSelection<String> accession = null;

	public DeleteExperimentTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		accession = new ListSingleSelection<>(new ArrayList<>(manager.getExperimentAccessions()));
	}

	public void run(TaskMonitor monitor) {
		experiment = accession.getSelectedValue();
		manager.deleteExperiment(experiment);
	}

	@ProvidesTitle
	public String title() { return "Delete Experiments"; }

	@Override
	public List<Class<?>> getResultClasses() {
		return Arrays.asList(String.class, JSONResult.class);
	}

	@Override
	public <R> R getResults(Class<? extends R> type) {
		if (type.equals(Experiment.class))
			return (R)experiment;
		else if (type.equals(JSONResult.class)) {
			JSONResult res = () -> {
				return "{}";
			};
			return (R)res;
		} else {
			return (R)("Removed experiment: "+experiment);
		}
	}
}

