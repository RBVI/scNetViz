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
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.MatrixUtils;

// Tunable to choose experiment?

public class HeatMapTask extends AbstractTask {
	final ScNVManager manager;
	final List<String> geneNames;
	final Category category;
	final Map<String, double[]> dataMap;
	final List<String> columnOrder;
	final boolean posOnly;
	int heatMapCount = -1;

	public HeatMapTask(final ScNVManager manager, final Category currentCategory, final List<String> rowLabels,
	                   final Map<String, double[]> dataMap, final List<String> columnOrder, boolean posOnly,
										 int count) {
		super();
		this.manager = manager;
		this.geneNames = rowLabels;
		this.category = currentCategory;
		this.dataMap = dataMap;
		this.columnOrder = columnOrder;
		this.posOnly = posOnly;
		this.heatMapCount = count;
	}

	// cyplot heat rowLabels="a,b,c,d,e,f" columnLabels="A,B,C" data="{\"A\":[1,2,3,4,5,6],\"B\":[-1,-2,-3,-4,-5,-6],\"C\":[0,1,-1,0,1,-1]}" title="Text Plot" xLabel="Upper" yLabel="Lower" editor=false
	public void run(TaskMonitor monitor) {
		Experiment exp = category.getExperiment();
		List<String> geneList = new ArrayList<>();
		String columnLabels = null;

		if (heatMapCount < 0)
			heatMapCount = Integer.parseInt(manager.getSetting(SETTING.HEATMAP_COUNT));

		for (String column: columnOrder) {
			double[] fc = dataMap.get(column);
			if (fc == null) continue;

			if (columnLabels == null)
				columnLabels = column;
			else
				columnLabels += ","+column;

			Integer[] sort = MatrixUtils.indexSort(fc, fc.length);

			// Skip over the NaN's
			int start = 0;
			for (start = 0; start < fc.length; start++) {
				if (!Double.isNaN(fc[sort[start]]))
					break;
			}

			int count = heatMapCount;
			if (!posOnly && heatMapCount < fc.length)
				count = heatMapCount/2;
			double[] topFC = new double[count];

			// System.out.println("count = "+count);
			// System.out.println("fc.length = "+fc.length);

			List<String> newGeneList = new ArrayList<String>();
			for (int topGene = 0; topGene < count; topGene++) {
				// System.out.println("Adding gene: "+geneNames.get(sort[topGene+start]));
				geneList.add(geneNames.get(sort[topGene+start]));
			}

			if (!posOnly && count < fc.length) {
				for (int topGene = fc.length-count; topGene < fc.length; topGene++) {
					// System.out.println("Adding gene: "+geneNames.get(sort[topGene]));
					geneList.add(geneNames.get(sort[topGene]));
				}
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
				// System.out.println(geneList.get(index)+" = "+fcData[index]);
			}
			sortedData.put(cat.toString(), fcData);
		}

		String data = CyPlotUtils.mapToData(sortedData);
		String accession = exp.getMetadata().get(Metadata.ACCESSION).toString();
		String title = exp.getSource().toString()+" "+ accession+ " Differential Expression";
		CyPlotUtils.createHeatMap(manager, CyPlotUtils.listToCSV(geneList), columnLabels, data, title, "Category", "Log(FC)", accession);
	}

}
