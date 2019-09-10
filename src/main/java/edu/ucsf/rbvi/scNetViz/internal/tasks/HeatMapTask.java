package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
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
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

// Tunable to choose experiment?

public class HeatMapTask extends AbstractTask {
	final ScNVManager manager;
	final Category category;
	final boolean posOnly;
	int heatMapCount = -1;
	final String selectedColumn;
	final DifferentialExpression diffExp;
	final double fdrCutoff;
	final double fcCutoff;
	final List<String> geneNames;
	final List<String> columnOrder;
	final Map<String, double[]> dataMap;

	public HeatMapTask(final ScNVManager manager, final Category currentCategory, 
	                   final DifferentialExpression diffExp, double fdrCutoff, double fcCutoff,
	                   boolean posOnly, int count,
	                   String selectedColumn) {

		super();
		this.manager = manager;
		this.posOnly = posOnly;
		this.heatMapCount = count;
		this.category = currentCategory;
		this.selectedColumn = selectedColumn;
		this.diffExp = diffExp;
		this.fdrCutoff = fdrCutoff;
		this.fcCutoff = fcCutoff;

		this.geneNames = diffExp.getRowLabels();
		this.dataMap = null;
		this.columnOrder = null;
	}

	public HeatMapTask(final ScNVManager manager, final Category currentCategory, final List<String> rowLabels,
	                   final Map<String, double[]> dataMap, final List<String> columnOrder, boolean posOnly,
	                   int count, final String column) {
		super();
		this.manager = manager;
		this.geneNames = rowLabels;
		this.category = currentCategory;
		this.dataMap = dataMap;
		this.columnOrder = columnOrder;
		this.posOnly = posOnly;
		this.heatMapCount = count;
		this.selectedColumn = column;

		this.diffExp = null;
		this.fdrCutoff = 1.00;
		this.fcCutoff = 0.0;
}


	// cyplot heat rowLabels="a,b,c,d,e,f" columnLabels="A,B,C" data="{\"A\":[1,2,3,4,5,6],\"B\":[-1,-2,-3,-4,-5,-6],\"C\":[0,1,-1,0,1,-1]}" title="Text Plot" xLabel="Upper" yLabel="Lower" editor=false
	public void run(TaskMonitor monitor) {
		Experiment exp = category.getExperiment();
		if (heatMapCount < 0)
			heatMapCount = Integer.parseInt(manager.getSetting(SETTING.HEATMAP_COUNT));

		List<String> geneList = new ArrayList<>();
		Map<String, double[]> dataMap;

		if (diffExp != null) {
			dataMap = getHeatMapFromDE(exp, geneList, monitor);
		} else {
			dataMap = getHeatMapFromMap(exp, geneList, monitor);
		}

		// Get our column labels
		String columnLabels = null;
		for (Object cat: diffExp.getLogGERMap().keySet()) {
			String column = category.mkLabel(cat);
			if (columnLabels == null)
				columnLabels = column;
			else
				columnLabels += ","+column;
		}

		String data = CyPlotUtils.mapToData(dataMap);
		String accession = exp.getMetadata().get(Metadata.ACCESSION).toString();
		String title = exp.getSource().toString()+" "+ accession+ " Differential Expression";
		CyPlotUtils.createHeatMap(manager, CyPlotUtils.listToCSV(geneList), columnLabels, 
		                          data, title, "Category", "Log(FC)", accession, posOnly);
	}

	private Map<String, double[]>  getHeatMapFromMap(Experiment exp, List<String> geneList, TaskMonitor monitor) {
		for (String column: columnOrder) {
			double[] fc = dataMap.get(column);
			if (fc == null) continue;

			if (selectedColumn != null && !column.equals(selectedColumn))
				continue;

			Integer[] sort = MatrixUtils.indexSort(fc, fc.length);

			// FIXME:  This really messes things up when we're getting called
			// with a small gene list -- i.e. when heatMapCount ~= fc.length
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

			List<String> newGeneList = new ArrayList<String>();
			for (int topGene = (fc.length-1); (topGene > (fc.length-count-1)) && (topGene >= 0); topGene--) {
				// System.out.println("Adding + gene: "+geneNames.get(sort[topGene])+"="+fc[sort[topGene]]);
				geneList.add(geneNames.get(sort[topGene]));
			}

			if (!posOnly && count < fc.length) {
				int first = start+count;
				if (first > fc.length) first = fc.length;
				for (int topGene = first; topGene > start; topGene--) {
					// System.out.println("Adding - gene: "+geneNames.get(sort[topGene])+"="+fc[sort[topGene]]);
					geneList.add(geneNames.get(sort[topGene+start]));
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
		return sortedData;
	}


	private Map<String, double[]> getHeatMapFromDE(Experiment exp, List<String> geneList, TaskMonitor monitor) {
		for (Object cat: diffExp.getLogGERMap().keySet()) {
			double[] fc = diffExp.getLogGER(cat, posOnly);
			if (fc == null) continue;

			String column = category.mkLabel(cat);

			if (selectedColumn != null && !column.equals(selectedColumn))
				continue;

			List<String> filteredGenes = new ArrayList<>();
			double[] filteredFC = diffExp.getGeneList(cat, fdrCutoff, fcCutoff, 0, posOnly, fc.length, filteredGenes);

			Integer[] sort = MatrixUtils.indexSort(filteredFC, filteredFC.length);

			/*
			// FIXME:  This really messes things up when we're getting called
			// with a small gene list -- i.e. when heatMapCount ~= fc.length
			// Skip over the NaN's
			int start = 0;
			for (start = 0; start < fc.length; start++) {
				if (!Double.isNaN(fc[sort[start]]))
					break;
			}
			*/

			// System.out.println("start = "+start);
			// System.out.println("fc.length = "+fc.length);
			// System.out.println("heatMapCount = "+heatMapCount);

			int count = heatMapCount;
			if (!posOnly && heatMapCount < filteredFC.length)
				count = heatMapCount/2;


			// System.out.println("count = "+count);
			// System.out.println("fc.length = "+fc.length);

			for (int topGene = (filteredFC.length-1); (topGene > (filteredFC.length-count-1)) && (topGene >= 0); topGene--) {
				// System.out.println("Adding + gene: "+filteredGenes.get(sort[topGene])+"="+filteredFC[sort[topGene]]);
				geneList.add(filteredGenes.get(sort[topGene]));
			}

			if (!posOnly && count < filteredFC.length) {
				int first = count;
				if (first > filteredFC.length) first = filteredFC.length;
				for (int topGene = first-1; topGene >= 0; topGene--) {
					// System.out.println("Adding - gene: "+filteredGenes.get(sort[topGene])+"="+filteredFC[sort[topGene]]);
					geneList.add(filteredGenes.get(sort[topGene]));
				}
			}
		}

		// Reverse the geneList
		Collections.reverse(geneList);

		// Now we have the list of genes that we want to use 
		Map<String, double[]> sortedData = new HashMap<>();
		for (Object cat: diffExp.getLogGERMap().keySet()) {
			double[] logGER = diffExp.getLogGER(cat, false);
			double[] fcData = new double[geneList.size()];
			for (int index = 0; index < geneList.size(); index++) {
				fcData[index] = logGER[geneNames.indexOf(geneList.get(index))];
			}
			sortedData.put(category.mkLabel(cat), fcData);
		}
		return sortedData;
	}

}
