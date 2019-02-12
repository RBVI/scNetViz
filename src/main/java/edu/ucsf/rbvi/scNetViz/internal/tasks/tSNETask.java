package edu.ucsf.rbvi.scNetViz.internal.tasks;

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

import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.BHTSne;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.FastTSne;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.MatrixOps;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.ParallelBHTsne;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.TSne;
import edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE.tSNEContext;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;


public class tSNETask extends AbstractTask implements ObservableTask {
	public static String SHORTNAME = "tsne";
	public static String NAME = "t-Distributed Stochastic Neighbor";
	public final static String GROUP_ATTRIBUTE = "__tSNE.SUID";

	@Tunable (description="Experiment to show")
	public ListSingleSelection<String> accession = null;

	@ContainsTunables
	public tSNEContext context = null;
	private DoubleMatrix matrix; 
	private double[][] tsneResult;
	final ScNVManager manager;
		
	public tSNETask(final ScNVManager manager) {
		this.manager = manager;
		this.context = new tSNEContext();
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
		matrix = null;
	}

	public tSNETask(DoubleMatrix matrix) {
		this.context = new tSNEContext();
		this.matrix = matrix;
		manager = null;
	}

	public String getShortName() { return SHORTNAME; }

	@ProvidesTitle
	public String getName() { return NAME; }
	
	public void run(TaskMonitor monitor) {
		monitor.setTitle(NAME);
		monitor.setStatusMessage("Running " + NAME);
		context.cancelled = false;
		MatrixOps matrixOps = new MatrixOps();
		
		if (matrix == null && accession == null) {
			monitor.showMessage(TaskMonitor.Level.ERROR, "Matrix is null");
			return;
		} else if (matrix == null) {
			Experiment experiment = manager.getExperiment(accession.getSelectedValue());
			matrix = (DoubleMatrix)experiment.getMatrix();
		}

		context.cancelled = false;

		// System.out.println("Getting matrix");
		List<String> rowLabels = new ArrayList<String>(matrix.getRowLabels());
		List<String> colLabels = new ArrayList<String>(matrix.getColLabels());
		// get sparse primitive matrix instead?
		double[][] Xin = matrix.getDoubleMatrix(0.0, false, true);
		// MatrixOps.debug("/tmp/originalMatrix", MatrixOps.doubleArrayToPrintString(matrix.getRowLabels(), Xin, false));

		long start = System.currentTimeMillis();
		
		Xin = MatrixOps.reduceMatrix(Xin, rowLabels, colLabels, 1, 1);
		// System.out.println("New size = "+Xin.length+"X"+Xin[0].length);
		// MatrixOps.debug("/tmp/reducedMatrix", MatrixOps.doubleArrayToPrintString(rowLabels, Xin, false));

		// Scale the data if we're supposed to
		if (context.logNormalize()) {
			// System.out.println("Normalizing matrix");
			Xin = MatrixOps.logNormalize(Xin, 10000);
			// MatrixOps.debug("/tmp/normalizedMatrix", MatrixOps.doubleArrayToPrintString(rowLabels, Xin, false));
		}

		if (context.centerAndScale()) {
			// System.out.println("Scaling matrix");
			Xin = MatrixOps.centerAndScale(Xin);
			// MatrixOps.debug("/tmp/scaledMatrix", MatrixOps.doubleArrayToPrintString(rowLabels, Xin, false));
		}

		// if (context.findVariableGenes()) {
		// System.out.println("Getting variable genes");
		BitSet variableGenes = new BitSet(Xin.length);
		Xin = MatrixOps.findVariableGenes(Xin, 0.0125, 3, 0.5, 20, variableGenes, rowLabels);
		// System.out.println("Xin["+Xin.length+"]["+Xin[0].length+"]");
		/*
		for (int row = 0; row < Xin.length; row++) {
			if (variableGenes.get(row)) {
				System.out.println("Variable gene = "+rowLabels.get(row));
			}
		}
		*/
		// }

		// Get the transposed matrix
		// Write the transposed matrix out
		Xin = matrixOps.transpose(Xin);
		// MatrixOps.debug("/tmp/transposedMatrix", MatrixOps.doubleArrayToPrintString(colLabels, Xin, false));
		context.setXin(Xin);

		context.setRowLabels(rowLabels);
		context.setColumnLabels(colLabels);

		monitor.showMessage(TaskMonitor.Level.INFO, "Setup complete in "+(System.currentTimeMillis()-start)/1000+" seconds");

		TSne tsne;
		try {
			if (context.useBarnesHut) {
				monitor.setTitle("Running t-Distributed Stochastic Neighbor (tSNE) using Barnes-Hut approximation");
				tsne = new ParallelBHTsne();
			} else {
				monitor.setTitle("Running t-Distributed Stochastic Neighbor (tSNE)");
				tsne = new FastTSne();
			}
		} catch (Exception e) {
			e.printStackTrace();
			monitor.showMessage(TaskMonitor.Level.ERROR, "Error calculating tSNE: "+e.getMessage());
			return;
		}

		tsneResult = tsne.tsne(context, monitor);
		if (tsneResult == null && context.cancelled) {
			monitor.setStatusMessage("Cancelled by user");
			return;
		}
		monitor.showMessage(TaskMonitor.Level.INFO, "tSNE complete in "+(System.currentTimeMillis()-start)/1000+" seconds");
	}

	public void cancel() {
		context.cancelled = true;
	}

	public double[][] getResults() { return tsneResult; }

	@Override
	public <R> R getResults(Class<? extends R> type) {
		if(tsneResult == null)
			return null;
		return (R) tsneResult;
	}
}
