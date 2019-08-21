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


// TODO: Consider random subsampling?
// TODO: Expose filter criteria
public class RemoteUMAPTask extends AbstractEmbeddingTask implements ObservableTask {
	public static String SHORTNAME = "umap";
	public static String NAME = "UMAP (Uniform Manifold Approximation and Projection)";
	public final static String GROUP_ATTRIBUTE = "__UMAP.SUID";

	@Tunable (description="Experiment to show")
	public ListSingleSelection<String> accession = null;

	@Tunable (description="Number of neighbors", 
	          longDescription="This parameter controls how UMAP balances local versus "+
	                          "global structure in the data."+
	                          "Low values will force UMAP to concentrate on very local structure, "+
	                          "while larger values will lose fine detail for the sake of getting "+
	                          "the broader context of the data.",
	          tooltip="<html>This parameter controls how UMAP balances local versus "+
	                   "global structure in the data."+
	                   "<p>Low values will force UMAP to concentrate on very local structure, "+
	                   "while larger values will lose fine detail for the sake of getting "+
	                   "the broader context of the data.</html>")
	public int n_neighbors = 10;

	@Tunable (description="Minimum distance",
	          tooltip="Controls how tightly UMAP is allowed to pack points together.  Lower numbers"+
	                          " result in tighter groupings and more distance between groups",
	          longDescription="Controls how tightly UMAP is allowed to pack points together. Lower numbers"+
	                          " result in tighter groupings and more distance between groups")
	public double min_dist = 0.5;

	public RemoteUMAPTask(final ScNVManager manager, final String acc) {
		super(manager);
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
		accession.setSelectedValue(acc);
	}

	public RemoteUMAPTask(final ScNVManager manager) {
		super(manager);
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
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
		String url = HTTPUtils.getWebServicesURL("umap", exp, "n_neighbors="+n_neighbors+"&min_dist="+min_dist);

		// Do the query
		try {
			List<String> lines = HTTPUtils.postFile(url, expFile, monitor);
			if (lines == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: UMAP failed: unable to read return");
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
			exp.setTSNE(embedding);
			exp.setPlotType("UMAP");
		} catch (Exception e) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: UMAP failed"+e.toString());
			return;
		}
		monitor.showMessage(TaskMonitor.Level.INFO, "UMAP complete in "+(System.currentTimeMillis()-start)/1000+" seconds");
	}
}
