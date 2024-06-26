package edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks;

import java.io.File;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileCategory;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileExperiment;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.view.CategoriesTab;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;

public class FileCategoryTask extends AbstractTask implements ObservableTask {
	final ScNVManager scManager;
	final FileSource fileSource;
	FileCategory cat = null;
	Experiment exp;

	@Tunable (description="Experiment to add category data to")
	public ListSingleSelection<Experiment> experiment = null;

	@Tunable (description="CSV file with category data",params="input=true")
	public File file;

	@Tunable (description="Data type of categories")
	public ListSingleSelection<String> dataType;

	@Tunable (description="Number of header columns (or rows if pivoted)")
	public int hdrCols=1;

	@Tunable (description="File needs to be pivoted (columns are categories)")
	public boolean pivot=false;

	@Tunable (description="Key column (or row, if pivoted)", format="#")
	public int keyColumn=0;

	public FileCategoryTask(final ScNVManager scManager, final FileSource fileSource, final Experiment exp) {
		super();
		this.scManager = scManager;
		this.fileSource = fileSource;
		this.exp = exp;

		if (exp != null)
			experiment = new ListSingleSelection<Experiment>(Collections.singletonList(exp));
		else {
			experiment = new ListSingleSelection<Experiment>(scManager.getExperiments());
		}
		dataType = new ListSingleSelection<String>("text","integer","float");
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		if (experiment != null)
			exp = experiment.getSelectedValue();

		taskMonitor.setTitle(getTitle());
		taskMonitor.setStatusMessage("Reading category file");

		try {
			cat = FileCategory.fetchCategory(scManager, exp, file,
			                                 dataType.getSelectedValue(), null,
			                                 pivot, hdrCols, keyColumn, false, taskMonitor);
      if (cat == null) {
		    taskMonitor.setStatusMessage("Failed to read category file: "+file);
        return;
      }
			exp.addCategory(cat);
			ExperimentFrame expFrame = scManager.getExperimentFrame(exp);
			if (expFrame != null) {
				String accession = (String)exp.getMetadata().get(Metadata.ACCESSION);
				CategoriesTab catTab = new CategoriesTab(scManager, exp, expFrame);
      	expFrame.addCategoriesContent(accession+": Categories Tab", catTab);
				catTab.changeCategory(cat, -1);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		taskMonitor.setStatusMessage("Successfully read category file");

	}

	@Override
	public <R> R getResults(Class<? extends R> type) {
		if (type.equals(Category.class))
			return (R)cat;
		return null;
	}

	@Override
	public List<Class<?>> getResultClasses() {
		return Arrays.asList(Category.class);
	}

	@ProvidesTitle
	public String getTitle() {return "Read category file";}
}
