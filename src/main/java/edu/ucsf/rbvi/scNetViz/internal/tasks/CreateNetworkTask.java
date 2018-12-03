package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

// Tunable to choose experiment?

public class CreateNetworkTask extends AbstractTask implements ObservableTask {
	final ScNVManager manager;

	// FIXME: these should be Tunables at some point
	double pValue;
	double log2FCCutoff;
	int nGenes;
	DifferentialExpression diffExp = null;

	public CreateNetworkTask(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public CreateNetworkTask(final ScNVManager manager, DifferentialExpression diffExp, 
	                         double pValue, double log2FCCutoff, int nGenes) {
		super();
		this.manager = manager;
		this.diffExp = diffExp;
		this.pValue = pValue;
		this.log2FCCutoff = log2FCCutoff;
		this.nGenes = nGenes;
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("Calculating Differential Expression");

		Category category = diffExp.getCurrentCategory();
		Set<Object> categoryValues = diffExp.getCategoryValues();

		List<String> allGenes = new ArrayList<String>();

		// Iterate over each category value
		for (Object cat: categoryValues) {
		//     Get the genes that match our criteria
		//     List<String> geneList = diffExp.getGeneList(cat, pValue, log2FCCutoff, nGenes);
		//     if (geneList != null && geneList.size() > 0) {
		//       allGenes.addAll(geneList);
		//       Create the network
		//       CyNetwork network = ModelUtils.createStringNetwork(manager, category.mkLabel(cat)+" Network", geneList);
		//       Create the columns
		//       ModelUtils.createDEColumns(manager, network, diffExp, category.mkLabel(cat));
		//       Add the data
		//       ModelUtils.updateDEData(manager, network, geneList, diffExp, category.mkLabel(cat));
		//       Style the network
		//       ModelUtils.addStyle(manager, network);
		//    }
		}
		//
		// Create the union network
		// Create the network
		// CyNetwork network = ModelUtils.createStringNetwork(manager, category.toString(), allGenes);
		for (Object cat: categoryValues) {
		//     Create the columns
		//     ModelUtils.createDEColumns(manager, network, diffExp, category.mkLabel(cat));
		//     Add the data
		//     ModelUtils.updateDEData(manager, network, allGenes, diffExp, category.mkLabel(cat));
		}
		// Style the network
		// ModelUtils.addStyle(manager, network);

	}

	// TODO: return the networks?
	public <R> R getResults(Class<? extends R> clazz) {
		if (clazz.equals(DifferentialExpression.class))
			return (R)diffExp;
		return null;
	}
}
