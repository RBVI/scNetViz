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
public class RemoteLeidenTask extends AbstractEmbeddingTask implements ObservableTask {
	public static String SHORTNAME = "leiden";
	public static String NAME = "Leiden clustering";
	public final static String GROUP_ATTRIBUTE = "__Leiden.SUID";
  String[] labels = {"Clusters"};

	@Tunable (description="Experiment to use")
	public ListSingleSelection<String> accession = null;

	@Tunable (description="Number of neighbors", 
	          longDescription="This parameter controls how Leiden balances local versus "+
	                          "global structure in the data."+
	                          "Low values will force Leiden to concentrate on very local structure, "+
	                          "while larger values will lose fine detail for the sake of getting "+
	                          "the broader context of the data.",
	          tooltip="<html>This parameter controls how Leiden balances local versus "+
	                   "global structure in the data."+
	                   "<p>Low values will force Leiden to concentrate on very local structure, "+
	                   "while larger values will lose fine detail for the sake of getting "+
	                   "the broader context of the data.</html>")
	public int n_neighbors = 15;

	@ContainsTunables
	public AdvancedRemoteParameters advancedParameters = null;

	public RemoteLeidenTask(final ScNVManager manager, final String acc) {
		super(manager);
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
		accession.setSelectedValue(acc);
		advancedParameters = new AdvancedRemoteParameters();
	}

	public RemoteLeidenTask(final ScNVManager manager) {
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

		Matrix mtx = exp.getMatrix();
		if ((mtx instanceof MatrixMarket)) {
      try {
        while (((MatrixMarket)mtx).cacheSent())
          Thread.sleep(1000);
      } catch (Exception e) {
        monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: Sleep interrupted!");
        return;
      }
    }

    /*
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
    */

		String url = HTTPUtils.getWebServicesURL("leiden", exp, 
		                                         "n_neighbors="+n_neighbors+
		                                         "&"+advancedParameters.getArgs());

		// Do the query
		try {
			List<String> lines = HTTPUtils.fetchResult(url, monitor);
			if (lines == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: Graph failed: unable to read return");
				return;
			}
			int lineNumber = 0;
			List<String[]> input = new ArrayList<>();
			String[] rowLabel = {"", "Leiden"};
			input.add(rowLabel);
			for (String line: lines) {
				String[] tokens = line.split(",");
				if (tokens.length < 2) continue;
				input.add(tokens);
				lineNumber++;
			}

			String categoryName = "Leiden Clusters";
			int index = 1;
			while (exp.getCategory(categoryName) != null) {
				categoryName = "Leiden Clusters - "+index;
				index++;
			}
			FileCategory louvainCategory = FileCategory.createCategory(manager, exp,
	                                                               categoryName, "integer", input,
	                                                               true, 0, 0, true, labels, monitor);
      try {
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
        e.printStackTrace();
        throw e;
      }
		} catch (Exception e) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: Calculation failed"+e.toString());
			return;
		}
		monitor.showMessage(TaskMonitor.Level.INFO, "Calculation complete in "+(System.currentTimeMillis()-start)/1000+" seconds");
	}
}
