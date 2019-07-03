package edu.ucsf.rbvi.scNetViz.internal.sources.hca.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCASource;

public class HCAShowEntriesTask extends AbstractTask {
	final ScNVManager scManager;
	final HCASource hcaSource;

	public HCAShowEntriesTask(final ScNVManager scManager, final HCASource hcaSource) {
		super();
		this.scManager = scManager;
		this.hcaSource = hcaSource;
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		taskMonitor.setTitle("Show Human Cell Atlas Entries Table");
		taskMonitor.setStatusMessage("Showing Human Cell Atlas Table");
		hcaSource.showEntriesTable(true);
	}

	@ProvidesTitle
	public String getTitle() {return "Show Human Cell Atlas Entries";}
}
