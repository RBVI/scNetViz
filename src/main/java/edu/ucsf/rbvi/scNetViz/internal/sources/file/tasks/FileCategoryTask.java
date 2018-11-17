package edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks;

import java.io.File;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileCategory;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileExperiment;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;

public class FileCategoryTask extends AbstractTask {
	final ScNVManager scManager;
	final FileSource fileSource;

	@Tunable (description="Experiment to add category data to")
	public ListSingleSelection<Experiment> experiment;

	@Tunable (description="CSV file with category data",params="input=true")
	public File file;

	@Tunable (description="Number of header columns (or rows if pivoted)")
	public int hdrCols=1;

	@Tunable (description="File needs to be pivoted (columns are categories)")
	public boolean pivot=false;

	public FileCategoryTask(final ScNVManager scManager, final FileSource fileSource) {
		super();
		this.scManager = scManager;
		this.fileSource = fileSource;
		experiment = new ListSingleSelection<Experiment>(scManager.getExperiments());
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		Experiment exp = experiment.getSelectedValue();

		taskMonitor.setTitle(getTitle());
		taskMonitor.setStatusMessage("Reading category file");

		try {
			FileCategory cat = FileCategory.fetchCategory(scManager, exp,
			                                              file, pivot, hdrCols, taskMonitor);
			exp.addCategory(cat);
		} catch (Exception e) {
			e.printStackTrace();
		}
		taskMonitor.setStatusMessage("Successfully read category file");

	}

	@ProvidesTitle
	public String getTitle() {return "Read category file";}
}
