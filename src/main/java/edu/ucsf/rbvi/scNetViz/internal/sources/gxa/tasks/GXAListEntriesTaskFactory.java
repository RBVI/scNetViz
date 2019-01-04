package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;

public class GXAListEntriesTaskFactory extends AbstractTaskFactory {
	final ScNVManager scManager;
	final GXASource gxaSource;

	public GXAListEntriesTaskFactory(final ScNVManager scManager, final GXASource gxaSource) {
		super();
		this.scManager = scManager;
		this.gxaSource = gxaSource;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator ti = new TaskIterator();
		if (gxaSource.getMetadata().size() == 0) {
			ti.append(new GXAFetchEntriesTask(scManager, gxaSource));
		}
		ti.append(new GXAListEntriesTask(scManager, gxaSource));
		return ti;
	}

	@Override
	public boolean isReady() { return true; }

}
