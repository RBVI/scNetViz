package edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;

public class FileCategoryTaskFactory extends AbstractTaskFactory {
	final ScNVManager scManager;
	final FileSource fileSource;
	final Experiment exp;

	public FileCategoryTaskFactory(final ScNVManager scManager, final FileSource fileSource) {
		super();
		this.scManager = scManager;
		this.fileSource = fileSource;
		this.exp = null;
	}

	public FileCategoryTaskFactory(final ScNVManager scManager, final FileSource fileSource,
	                               final Experiment exp) {
		super();
		this.scManager = scManager;
		this.fileSource = fileSource;
		this.exp = exp;
	}

	public TaskIterator createTaskIterator() {
		return new TaskIterator(new FileCategoryTask(scManager, fileSource, exp));
	}

	public boolean isReady() {
		if (scManager.getExperiments() == null || scManager.getExperiments().size() == 0) 
			return false;
		return true;
	}
}
