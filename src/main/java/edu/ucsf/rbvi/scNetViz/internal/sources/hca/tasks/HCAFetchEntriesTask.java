package edu.ucsf.rbvi.scNetViz.internal.sources.hca.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCASource;

public class HCAFetchEntriesTask extends AbstractTask {
	final ScNVManager scManager;
	final HCASource hcaSource;

	public HCAFetchEntriesTask(final ScNVManager scManager, final HCASource hcaSource) {
		super();
		this.scManager = scManager;
		this.hcaSource = hcaSource;
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		taskMonitor.setTitle(getTitle());
		taskMonitor.setStatusMessage("Fetching all HCA entries");
		hcaSource.loadHCAEntries(taskMonitor);
	}

	@ProvidesTitle
	public String getTitle() {return "HCA Fetch Entries";}
}
