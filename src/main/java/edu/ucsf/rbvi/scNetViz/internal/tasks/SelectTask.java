package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;


public class SelectTask extends AbstractTask {

	@Tunable(description="Genes to select", context="nogui")
	public String genes = null;

	@Tunable(description="Cells to select", context="nogui")
	public String cells = null;

	@Tunable(description="Accession of experiment", context="nogui")
	public String accession = null;

	final ScNVManager manager;

	public SelectTask(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("Selecting genes and assays");

		List<Experiment> exps = manager.getExperiments();
		if (exps.size() == 1 && accession == null)
			accession = exps.get(0).getMetadata().get(Metadata.ACCESSION).toString();

		if (accession == null) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Experiment accession must be specified!");
			return;
		}

		Experiment exp = manager.getExperiment(accession);
		ExperimentFrame expFrame = manager.getExperimentFrame(exp);
		if (genes != null) {
			List<String> geneList = getListFromString(genes);
			expFrame.selectGenes(geneList);
		}

		if (cells != null) {
			List<String> assayList = getListFromString(cells);
			expFrame.selectAssays(assayList);
		}
	}

	List<String> getListFromString(String csvString) {
		if (csvString == null) return null;
		String[] splitString = csvString.split("\t");
		if (splitString.length == 1)
			splitString = csvString.split(",");

		List<String> arrayString = new ArrayList<>();
		for (String str: splitString) {
			arrayString.add(str);
		}
		return arrayString;
	}

}
