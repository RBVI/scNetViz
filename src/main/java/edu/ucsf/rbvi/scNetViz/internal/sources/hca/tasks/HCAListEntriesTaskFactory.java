package edu.ucsf.rbvi.scNetViz.internal.sources.hca.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCASource;

public class HCAListEntriesTaskFactory extends AbstractTaskFactory {
	final ScNVManager scManager;
	final HCASource hcaSource;

	public HCAListEntriesTaskFactory(final ScNVManager scManager, final HCASource hcaSource) {
		super();
		this.scManager = scManager;
		this.hcaSource = hcaSource;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator ti = new TaskIterator();
		if (hcaSource.getMetadata().size() == 0) {
			ti.append(new HCAFetchEntriesTask(scManager, hcaSource));
		}
		ti.append(new HCAListEntriesTask(scManager, hcaSource));
		return ti;
	}

	@Override
	public boolean isReady() { return true; }

}
