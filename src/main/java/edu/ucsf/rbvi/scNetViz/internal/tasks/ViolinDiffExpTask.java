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

public class ViolinDiffExpTask extends AbstractTask {
	final ScNVManager manager;
	final List<String> geneNames;
	final Category category;
	final Map<String, double[]> dataMap;
	final List<String> columnOrder;
	final boolean posOnly;
	final String title;
	final String accession;

	public ViolinDiffExpTask(final ScNVManager manager, final Category currentCategory, 
	                         final DifferentialExpression diffExp, boolean posOnly) {
		this.manager = manager;
		this.posOnly = posOnly;
		this.category = currentCategory;


		dataMap = new HashMap<>();
		columnOrder = new ArrayList<>();
		for (Object cat: diffExp.getLogGERMap().keySet()) {
			double[] logGER = diffExp.getLogGER(cat, posOnly);
			if (logGER != null) {
				dataMap.put(currentCategory.mkLabel(cat), logGER);
				columnOrder.add(currentCategory.mkLabel(cat));
			}
		}
		geneNames = diffExp.getRowLabels();

		accession = diffExp.getExperiment().getMetadata().get(Metadata.ACCESSION).toString();
		String source = diffExp.getExperiment().getSource().toString();
		title = source.toString()+" "+ accession+ " Differential Expression";
	}

	// cyplot heat rowLabels="a,b,c,d,e,f" columnLabels="A,B,C" data="{\"A\":[1,2,3,4,5,6],\"B\":[-1,-2,-3,-4,-5,-6],\"C\":[0,1,-1,0,1,-1]}" title="Text Plot" xLabel="Upper" yLabel="Lower" editor=false
	public void run(TaskMonitor monitor) {
		String[] dataAndNames = CyPlotUtils.mapToDataAndNames(dataMap, geneNames, columnOrder);
		CyPlotUtils.createViolinPlot(manager, dataAndNames[0], dataAndNames[1], 
		                             CyPlotUtils.listToCSV(columnOrder), title, "", "Log(FC)", accession, false, 0.7);
	}

}
