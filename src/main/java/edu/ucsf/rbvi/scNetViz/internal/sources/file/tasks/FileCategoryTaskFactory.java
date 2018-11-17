package edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;

public class FileCategoryTaskFactory extends AbstractTaskFactory {
	final ScNVManager scManager;
	final FileSource fileSource;

	public FileCategoryTaskFactory(final ScNVManager scManager, final FileSource fileSource) {
		super();
		this.scManager = scManager;
		this.fileSource = fileSource;
	}

	public TaskIterator createTaskIterator() {
		return new TaskIterator(new FileCategoryTask(scManager, fileSource));
	}

	public boolean isReady() {
		if (scManager.getExperiments() == null || scManager.getExperiments().size() == 0) 
			return false;
		return true;
	}
}
