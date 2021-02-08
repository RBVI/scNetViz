package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
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
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.MatrixUtils;

// Tunable to choose experiment?

public class ViolinGeneTask extends AbstractTask {
	final ScNVManager manager;
	final Category category;
	final int selectedCategory;
	final int selectedGene;
	final Map<String, double[]> dataMap;
	final List<String> columnOrder;
	final String title;
	final String accession;
	final Experiment experiment;
	final List<String> cellNames;

	public ViolinGeneTask(final ScNVManager manager, final Category currentCategory, 
	                      final int selectedCategory, final int selectedGene) {
		this.manager = manager;
		this.category = currentCategory;
		this.experiment = currentCategory.getExperiment();
		if (selectedCategory == -1) 
			this.selectedCategory = currentCategory.getDefaultRow();
		else
			this.selectedCategory = selectedCategory;

		this.selectedGene = selectedGene;

		Map<Object, List<Integer>> catMap = currentCategory.getCatMap(selectedCategory);

		dataMap = new HashMap<>();
		columnOrder = new ArrayList<>();
		cellNames = experiment.getMatrix().getColLabels(0);
		for (Object cat: catMap.keySet()) {
			// System.out.println("Cat: "+cat.toString());
			double[] tpm = new double[cellNames.size()];
			Arrays.fill(tpm, Double.NaN);
			for (int cell: catMap.get(cat)) {
				Object v = experiment.getTableModel().getValueAt(selectedGene, cell+1);
				if (v == null) 
					tpm[cell] = 0.0;
				else
					tpm[cell] = (Double)v;
				// System.out.println("selectedGene = "+selectedGene+", cell = "+(cell+1)+" = "+tpm[cell]);
			}
			dataMap.put(currentCategory.mkLabel(cat), tpm);
			columnOrder.add(currentCategory.mkLabel(cat));
		}
		accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
		String source = experiment.getSource().toString();
		title = source.toString()+" "+ accession+ " "+experiment.getMatrix().getRowLabel(selectedGene);
	}

	// cyplot heat rowLabels="a,b,c,d,e,f" columnLabels="A,B,C" data="{\"A\":[1,2,3,4,5,6],\"B\":[-1,-2,-3,-4,-5,-6],\"C\":[0,1,-1,0,1,-1]}" title="Text Plot" xLabel="Upper" yLabel="Lower" editor=false
	public void run(TaskMonitor monitor) {
		String[] dataAndNames = CyPlotUtils.mapToDataAndNames(dataMap, cellNames, columnOrder);
		CyPlotUtils.createViolinPlot(manager, dataAndNames[0], dataAndNames[1], 
		                             CyPlotUtils.listToCSV(columnOrder), title, "", "TPM", accession, false, 0.7);
	}

}
