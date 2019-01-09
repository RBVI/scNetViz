package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ContainsTunables;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.BHTSne;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.FastTSne;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.TSne;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.tSNEContext;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;


public class tSNETask extends AbstractTask implements ObservableTask {
	public static String SHORTNAME = "tsne";
	public static String NAME = "t-Distributed Stochastic Neighbor";
	public final static String GROUP_ATTRIBUTE = "__tSNE.SUID";

	@ContainsTunables
	public tSNEContext context = null;
	private DoubleMatrix matrix; 
	private double[][] tsneResult;
		
	public tSNETask(DoubleMatrix matrix) {
		this.context = new tSNEContext();
		this.matrix = matrix;
	}

	public String getShortName() { return SHORTNAME; }

	@ProvidesTitle
	public String getName() { return NAME; }
	
	public void run(TaskMonitor monitor) {
		monitor.setTitle(NAME);
		monitor.setStatusMessage("Running " + NAME);
		context.cancelled = false;
		
		if (matrix == null) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Matrix is null");
			return;
		}
		
		context.cancelled = false;
		context.setXin(matrix.getDoubleMatrix(0.0));
		
		TSne tsne;
		if (context.useBarnesHut) {
			monitor.setTitle("Running t-Distributed Stochastic Neighbor (tSNE) using Barnes-Hut approximation");
			tsne = new BHTSne();
		} else {
			monitor.setTitle("Running t-Distributed Stochastic Neighbor (tSNE)");
			tsne = new FastTSne();
		}

		tsneResult = tsne.tsne(context, monitor);
		if (tsneResult == null && context.cancelled) {
			monitor.setStatusMessage("Cancelled by user");
			return;
		}
	}

	public void cancel() {
		context.cancelled = true;
	}

	@Override
	public <R> R getResults(Class<? extends R> type) {
		if(tsneResult == null)
			return null;
		return (R) tsneResult;
	}
}
