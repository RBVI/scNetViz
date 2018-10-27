package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;

public class GXAShowEntriesTask extends AbstractTask {
	final ScNVManager scManager;
	final GXASource gxaSource;

	public GXAShowEntriesTask(final ScNVManager scManager, final GXASource gxaSource) {
		super();
		this.scManager = scManager;
		this.gxaSource = gxaSource;
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		taskMonitor.setTitle("Show GXA Entries Table");
		taskMonitor.setStatusMessage("Showing Entries Table");
		gxaSource.showEntriesTable(true);
	}

	@ProvidesTitle
	public String getTitle() {return "Show GXA Entries";}
}
