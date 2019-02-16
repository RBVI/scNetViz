package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks;

import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ShowExperimentTableTask;


public class GXALoadExperimentTask extends AbstractTask {
	final ScNVManager scManager;
	final GXASource gxaSource;

	@Tunable(description="Experiment accession id", required=true)
	public String accession;

	@Tunable (description="Show experiment table after loading",context="nogui")
	public boolean showTable = true;

	public GXALoadExperimentTask(final ScNVManager scManager, final GXASource gxaSource) {
		super();
		this.scManager = scManager;
		this.gxaSource = gxaSource;
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		if (accession == null) {
			taskMonitor.showMessage(TaskMonitor.Level.ERROR, "No accession id provided");
			return;
		}
		taskMonitor.setTitle("Loading Single Cell Expression Atlas Experiment: "+accession);
		List<Metadata> metadata = gxaSource.getMetadata();
		Experiment experiment = gxaSource.getExperiment(accession, taskMonitor, showTable);

		if (showTable) {
			insertTasksAfterCurrentTask(new ShowExperimentTableTask(scManager, experiment));
		}
	}

	@ProvidesTitle
	public String getTitle() {return "Load Single Cell Expression Atlas Experiment";}
}
