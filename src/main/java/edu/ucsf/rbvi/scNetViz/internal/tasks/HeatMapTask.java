package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.MatrixUtils;

// Tunable to choose experiment?

public class HeatMapTask extends AbstractTask {
	final ScNVManager manager;
	final DifferentialExpression diffExp;
	final Category category;
	final Map<String, double[]> dataMap;

	public HeatMapTask(final ScNVManager manager, DifferentialExpression diffExp, Category currentCategory,
	                   final Map<String, double[]> dataMap) {
		super();
		this.manager = manager;
		this.diffExp = diffExp;
		this.category = currentCategory;
		this.dataMap = dataMap;
	}

	// cyplot heat rowLabels="a,b,c,d,e,f" columnLabels="A,B,C" data="{\"A\":[1,2,3,4,5,6],\"B\":[-1,-2,-3,-4,-5,-6],\"C\":[0,1,-1,0,1,-1]}" title="Text Plot" xLabel="Upper" yLabel="Lower" editor=false
	public void run(TaskMonitor monitor) {
		Experiment exp = diffExp.getExperiment();
		List<String> geneNames = exp.getMatrix().getRowLabels();
		List<String> geneList = new ArrayList<>();
		String columnLabels = null;

		for (Object cat: dataMap.keySet()) {
			double[] fc = dataMap.get(cat);
			if (fc == null) continue;

			if (columnLabels == null)
				columnLabels = cat.toString();
			else
				columnLabels += ","+cat.toString();

			Integer[] sort = MatrixUtils.indexSort(fc, fc.length);

			// Skip over the NaN's
			int start = 0;
			for (start = 0; start < fc.length; start++) {
				if (!Double.isNaN(fc[sort[start]]))
					break;
			}

			double[] topFC = new double[20];
			// Now get the top 10 and the bottom 10
			List<String> newGeneList = new ArrayList<String>();
			for (int topGene = 0; topGene < 10; topGene++) {
				System.out.println("Top gene("+(topGene+start)+"): "+diffExp.getRowLabel(sort[topGene+start])+" = "+fc[sort[topGene+start]]);
				geneList.add(diffExp.getRowLabel(sort[topGene+start]));
			}

			for (int topGene = fc.length-10; topGene < fc.length-1; topGene++) {
				System.out.println("Bottom gene("+(topGene)+"): "+diffExp.getRowLabel(sort[topGene])+" = "+fc[sort[topGene]]);
				geneList.add(diffExp.getRowLabel(sort[topGene]));
			}
		}

		// Reverse the geneList
		Collections.reverse(geneList);

		// Now we have the list of genes that we want to use 
		Map<String, double[]> sortedData = new HashMap<>();
		for (Object cat: dataMap.keySet()) {
			double[] fc = dataMap.get(cat);
			double[] fcData = new double[geneList.size()];
			for (int index = 0; index < geneList.size(); index++) {
				fcData[index] = fc[geneNames.indexOf(geneList.get(index))];
			}
			sortedData.put(cat.toString(), fcData);
		}

		String data = CyPlotUtils.mapToData(sortedData);
		String accession = exp.getMetadata().get(Metadata.ACCESSION).toString();
		String title = exp.getSource().toString()+" "+ accession+ " Differential Expression";
		CyPlotUtils.createHeatMap(manager, CyPlotUtils.listToCSV(geneList), columnLabels, data, title, "Category", "Log(FC)", accession);
	}

}
