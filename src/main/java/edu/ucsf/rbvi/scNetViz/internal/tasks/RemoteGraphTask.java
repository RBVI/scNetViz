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
public class RemoteGraphTask extends AbstractEmbeddingTask implements ObservableTask {
	public static String SHORTNAME = "draw_graph";
	public static String NAME = "Force-directed graph drawing";
	public final static String GROUP_ATTRIBUTE = "__Graph.SUID";
	public String[] layouts = {"fa (ForceAtlas2)","kk (Kamada Kawai)", "fr (Fruchterman Reingold)",
	                           "lgl (Large Graph)","dlr (Distributed Recursive Layout)",
	                           "rt (Reingold Tilford tree layout)"};

	@Tunable (description="Experiment to show")
	public ListSingleSelection<String> accession = null;

	@Tunable (description="Graph layout algorithm", 
	          longDescription="‘fa’ (ForceAtlas2) or any valid igraph layout. Of particular interest are "+
	                          "‘fr’ (Fruchterman Reingold), ‘kk’ (Kamadi Kawai’, slower than ‘fr’), ‘lgl’ "+
	                          "(Large Graph, very fast), ‘drl’ (Distributed Recursive Layout, pretty fast) and "+
	                          "‘rt’ (Reingold Tilford tree layout).",
	          tooltip="<html><body>Available options are: "+
	                  "<ul><li>‘fa’ (ForceAtlas2)</li>"+
	                  "<li>‘fr’ (Fruchterman Reingold)</li>"+
	                  "<li>‘kk’ (Kamadi Kawai, slower than ‘fr’)</li>"+
	                  "<li>‘lgl’ (Large Graph, very fast)</li>"+
	                  "<li>‘drl’ (Distributed Recursive Layout, pretty fast)</li>"+
	                  "<li>‘rt’ (Reingold Tilford tree layout)</li></ul></body></html>")
	public ListSingleSelection<String> layout = new ListSingleSelection<>(layouts);

	@ContainsTunables
	public AdvancedRemoteParameters advancedParameters = null;

	public RemoteGraphTask(final ScNVManager manager, final String acc) {
		super(manager);
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
		accession.setSelectedValue(acc);
	}

	public RemoteGraphTask(final ScNVManager manager) {
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

		// Split off the explanatory text
		String type = layout.getSelectedValue().split(" ")[0];

		String url = HTTPUtils.getWebServicesURL("drawgraph", exp, "layout="+type+
		                                         "&"+advancedParameters.getArgs());

		// Do the query
		try {
			List<String> lines = HTTPUtils.postFile(url, expFile, monitor);
			if (lines == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: Graph failed: unable to read return");
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
			exp.setPlotType("Graph");
		} catch (Exception e) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "ERROR: Graph failed"+e.toString());
			return;
		}
		monitor.showMessage(TaskMonitor.Level.INFO, "Graph complete in "+(System.currentTimeMillis()-start)/1000+" seconds");
	}
}
