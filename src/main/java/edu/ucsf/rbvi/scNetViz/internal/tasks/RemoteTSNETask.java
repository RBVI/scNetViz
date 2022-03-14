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
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;


public class RemoteTSNETask extends AbstractEmbeddingTask {
	public static String SHORTNAME = "tsne";
	public static String NAME = "t-Distributed Stochastic Neighbor";
	public final static String GROUP_ATTRIBUTE = "__TSNE.SUID";

	@Tunable (description="Experiment to show")
	public ListSingleSelection<String> accession = null;

	@Tunable (description="Perplexity", 
	          tooltip="<html>Related to the number of nearest neighbors.  Larger datasets usually require"+
	                  "<br/>a larger perplexity. Consider selecting a value between 5 and 50.</html>",
	          longDescription="Related to the number of nearest neighbors.  Larger datasets usually require"+
	                          " a larger perplexity. Consider selecting a value between 5 and 50.")
	public double perplexity = 20;

	@Tunable (description="Initial dimensions", 
	          tooltip="The number of principal components to use.  -1 chooses is the same as scanpy's"+
	                  "	'None', and 0 means don't do PCA.",
	          format="#",
	          longDescription="The number of principal components to use.")
	public int n_pcs = -1;

	@Tunable (description="Early exaggeration", 
	          tooltip="<html>Controls how tight natural clusters in the original space "+
	                  "are in the embedded space and how much space will be between "+
	                  "them. <br/>For larger values, the space between natural clusters will be "+
	                  "larger in the embedded space.</html>",
	          longDescription="Controls how tight natural clusters in the original space "+
	                          "are in the embedded space and how much space will be between "+
	                          "them. For larger values, the space between natural clusters will be "+
	                          "larger in the embedded space.") 
	public double early_exaggeration = 12.0;

	@Tunable (description="Learning rate", 
	          tooltip="<html>The learning rate can be a critical parameter. It should <br/>"+
	                  "be between 100 and 1000. If the cost function increases during <br/>"+
	                  "initial optimization, the early exaggeration factor or the <br/>"+
	                  "learning rate might be too high. If the cost function gets <br/>"+
	                  "stuck in a bad local minimum increasing the learning rate helps sometimes.</html>",
	          longDescription="The learning rate can be a critical parameter. It should "+
	                          "be between 100 and 1000. If the cost function increases during "+
	                          "initial optimization, the early exaggeration factor or the "+
	                          "learning rate might be too high. If the cost function gets "+
	                          "stuck in a bad local minimum increasing the learning rate helps sometimes.")
	public double learning_rate = 1000.0;

	@ContainsTunables
	public AdvancedRemoteParameters advancedParameters = null;

	public RemoteTSNETask(final ScNVManager manager, final String acc) {
		super(manager);
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
		accession.setSelectedValue(acc);
		advancedParameters = new AdvancedRemoteParameters();
	}

	public RemoteTSNETask(final ScNVManager manager) {
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

		String query = "perplexity="+perplexity+"&early_exaggeration="+early_exaggeration;
		if (n_pcs >= 0) query += "&n_pcs="+n_pcs;
		query += "&"+advancedParameters.getArgs();
		String url = HTTPUtils.getWebServicesURL("tsne", exp, query);

		// Do the query
		try {
			// List<String> lines = HTTPUtils.postFile(url, expFile, monitor);
			List<String> lines = HTTPUtils.fetchResult(url, monitor);
			if (lines == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: TSNE failed: unable to read return");
				return;
			}
			embedding = new double[lines.size()-1][2];
			int lineNumber = 0;
			for (String line: lines) {
				String[] tokens = line.split(",");
				if (tokens.length <= 2) continue;
				embedding[lineNumber][0] = Double.valueOf(tokens[tokens.length-2]);
				embedding[lineNumber][1] = Double.valueOf(tokens[tokens.length-1]);
				lineNumber++;
			}
			scale(embedding); // Scale everything between 0 and 1 so that it appears unitless
      updateCellPlot(exp, "tSNE");
		} catch (Exception e) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: t-SME failed"+e.toString());
			return;
		}
		monitor.showMessage(TaskMonitor.Level.INFO, "t-SME complete in "+(System.currentTimeMillis()-start)/1000+" seconds");
	}

}
