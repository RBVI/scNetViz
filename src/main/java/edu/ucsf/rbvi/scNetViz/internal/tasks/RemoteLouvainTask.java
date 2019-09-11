package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ContainsTunables;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileCategory;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.CategoriesTab;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;


// TODO: Consider random subsampling?
// TODO: Expose filter criteria
public class RemoteLouvainTask extends AbstractEmbeddingTask implements ObservableTask {
	public static String SHORTNAME = "louvain";
	public static String NAME = "Louvain clustering";
	public final static String GROUP_ATTRIBUTE = "__Louvain.SUID";

	@Tunable (description="Experiment to use")
	public ListSingleSelection<String> accession = null;

	@Tunable (description="Number of neighbors", 
	          longDescription="This parameter controls how Louvain balances local versus "+
	                          "global structure in the data."+
	                          "Low values will force Louvain to concentrate on very local structure, "+
	                          "while larger values will lose fine detail for the sake of getting "+
	                          "the broader context of the data.",
	          tooltip="<html>This parameter controls how Louvain balances local versus "+
	                   "global structure in the data."+
	                   "<p>Low values will force Louvain to concentrate on very local structure, "+
	                   "while larger values will lose fine detail for the sake of getting "+
	                   "the broader context of the data.</html>")
	public int n_neighbors = 15;

	@ContainsTunables
	public AdvancedRemoteParameters advancedParameters = null;

	public RemoteLouvainTask(final ScNVManager manager, final String acc) {
		super(manager);
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
		accession.setSelectedValue(acc);
		advancedParameters = new AdvancedRemoteParameters();
	}

	public RemoteLouvainTask(final ScNVManager manager) {
		super(manager);
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
		advancedParameters = new AdvancedRemoteParameters();
	}

	public String getShortName() { return SHORTNAME; }

	@ProvidesTitle
	public String getName() { return NAME; }
	
	public void run(TaskMonitor monitor) {
		monitor.setTitle(NAME);
		monitor.setStatusMessage("Running " + NAME + " on server");
		long start = System.currentTimeMillis();

		// Get the experiment
		Experiment exp = manager.getExperiment(accession.getSelectedValue());
		// Get the MatrixMarket matrix
		Matrix mtx = exp.getMatrix();
		if (!(mtx instanceof MatrixMarket)) {
			monitor.showMessage(TaskMonitor.Level.ERROR,"Matrix must be of type MatrixMarket");
			return;
		}
		MatrixMarket mmtx = (MatrixMarket)mtx;

		if (!mmtx.hasCache()) {
			createCache(mmtx, exp);
		}
		File expFile = mmtx.getMatrixCache();

		String url = HTTPUtils.getWebServicesURL("louvain", exp, 
		                                         "n_neighbors="+n_neighbors+
		                                         "&"+advancedParameters.getArgs());

		// Do the query
		try {
			List<String> lines = HTTPUtils.postFile(url, expFile, monitor);
			if (lines == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: Graph failed: unable to read return");
				return;
			}
			int lineNumber = 0;
			List<String[]> input = new ArrayList<>(lines.size());
			String[] rowLabel = {"", "Louvain"};
			input.add(rowLabel);
			for (String line: lines) {
				String[] tokens = line.split(",");
				if (tokens.length < 2) continue;
				input.add(tokens);
				lineNumber++;
			}
			String categoryName = "Louvain Clusters";
			int index = 1;
			while (exp.getCategory(categoryName) != null) {
				categoryName = "Louvain Clusters - "+index;
				index++;
			}
			FileCategory louvainCategory = FileCategory.createCategory(manager, exp,
	                                                               categoryName, "integer", input,
	                                                               true, 1, true, monitor);
			exp.addCategory(louvainCategory);
			ExperimentFrame expFrame = manager.getExperimentFrame(exp);
			if (expFrame != null) {
				String accession = (String)exp.getMetadata().get(Metadata.ACCESSION);
				CategoriesTab catTab = expFrame.getCategoriesTab();
				catTab = new CategoriesTab(manager, exp, expFrame);
      	expFrame.addCategoriesContent(accession+": Categories Tab", catTab);
				catTab.changeCategory(louvainCategory, -1);
			}
		} catch (Exception e) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: Calculation failed"+e.toString());
			return;
		}
		monitor.showMessage(TaskMonitor.Level.INFO, "Calculation complete in "+(System.currentTimeMillis()-start)/1000+" seconds");
	}
}