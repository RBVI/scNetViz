package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;

public class GXAFetchEntriesTask extends AbstractTask {
	final ScNVManager scManager;
	final GXASource gxaSource;

	public GXAFetchEntriesTask(final ScNVManager scManager, final GXASource gxaSource) {
		super();
		this.scManager = scManager;
		this.gxaSource = gxaSource;
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		taskMonitor.setTitle(getTitle());
		taskMonitor.setStatusMessage("Fetching all GXA entries");
		gxaSource.loadGXAEntries(taskMonitor);
	}

	@ProvidesTitle
	public String getTitle() {return "GXA Fetch Entries";}
}
