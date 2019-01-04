package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.json.JSONResult;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;

public class ListExperimentsTask extends AbstractTask implements ObservableTask { 
	final ScNVManager manager;

	public ListExperimentsTask(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("List experiments");
	}

	@Override
	public List<Class<?>> getResultClasses() {
		return Arrays.asList(String.class, JSONResult.class, List.class);
	}

	@Override
	public <R> R getResults(Class<? extends R> type) {
		List<Experiment> experiments = manager.getExperiments();
		if (type.equals(List.class)) {
			return (R) experiments;
		} else if (type.equals(JSONResult.class)) {
			return (R) (new JSONExp(experiments));
		} else {
			String str = "Current experiments: \n";
			for (Experiment exp: experiments) {
				str += "    "+exp.getMetadata().get(Metadata.ACCESSION).toString()+": ";
				str += exp.getMetadata().get(Metadata.DESCRIPTION).toString();
				str += " ["+exp.getMatrix().getNRows()+"x"+exp.getMatrix().getNCols()+"]\n";
			}
			return (R) str;
		}
	}

	class JSONExp implements JSONResult {
		final List<Experiment> experiments;

		JSONExp(final List<Experiment> exps) {
			experiments = exps;
		}

		public String getJSON() {
			StringBuilder builder = new StringBuilder();
			builder.append("[");
			for (Experiment exp: experiments) {
				builder.append("{");
				builder.append("source: '"+exp.getSource().toString()+"',");
				builder.append("accession: '"+exp.getMetadata().get(Metadata.ACCESSION).toString()+"',");
				builder.append("species: '"+exp.getSpecies().toString()+"',");
				builder.append("description: '"+exp.getMetadata().get(Metadata.DESCRIPTION).toString()+"',");
				builder.append("rows: '"+exp.getMatrix().getNRows()+"',");
				builder.append("columns: '"+exp.getMatrix().getNCols()+"',");
				builder.append("categories: '"+exp.getCategories().size()+"'}");
			}
			return builder.substring(0, builder.length()-1)+"]";
		}
	}

}
