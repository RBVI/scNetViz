package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;
import java.util.Arrays;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;


public class SettingsTask extends AbstractTask {

	@Tunable(description="Log2FC cutoff", groups={"Differential Expression Calculation Settings"})
	public double deLogFCCutoff = 0.5;

	@Tunable(description="Min.pct cutoff", groups={"Differential Expression Calculation Settings"})
	public double deMinPctCutoff = 10.0;

	@Tunable(description="Log2FC cutoff", groups={"Network Creation Settings"})
	public double netLogFCCutoff = 1.0;

	@Tunable(description="pValue cutoff", groups={"Network Creation Settings"})
	public double netPvCutoff = 0.05;

	@Tunable(description="Top n genes", 
	         tooltip="<html>This will select the <i>n</i> most significant genes</html>",
	         groups={"Network Creation Settings"})
	public int topGenes = 0;

	@Tunable(description="Maximum Genes", 
	         tooltip="<html>Limits the number of genes per category to use for network creation</html>",
	         groups={"Network Creation Settings"})
	public int maxGenes = 50;

	@Tunable(description="Positive only", 
	         tooltip="<html>Only include genes with a positive Log2FC</html>",
	         groups={"Network Creation Settings"})
	public boolean positiveOnly = false;

	@Tunable(description="Heat map count limit", 
	         tooltip="<html>The number of top genes for each category will be limited by this amount.</html>",
	         groups={"Visualization Settings"})
	public int heatMapCount = 20;

	@Tunable(description="Double-Click Action", 
	         tooltip="<html>By default, don't automatically create the networks</html>")
	public ListSingleSelection<String> doubleClickAction = new ListSingleSelection<String>(Arrays.asList("View Data", "Create Network"));

	final ScNVManager manager;

	/*
		MAX_GENES("maxGenes", "500"),
		DE_FC_CUTOFF("deFcCutoff", "0.5"),
		DE_MIN_PCT_CUTOFF("deMinPctCutoff", "10"),
		TOP_GENES("topGenes",""),
		NET_PV_CUTOFF("netPvCutoff","0.05"),
		NET_FC_CUTOFF("netFcCutoff","1.00"),
		DONT_ANALYZE("dontAnalyze","false");
	*/

	public SettingsTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		maxGenes = Integer.parseInt(manager.getSetting(SETTING.MAX_GENES));
		String strTopGenes = manager.getSetting(SETTING.TOP_GENES);
		if (strTopGenes == null || strTopGenes.length() == 0)
			topGenes = 0;
		else
			topGenes = Integer.parseInt(manager.getSetting(SETTING.TOP_GENES));
		deLogFCCutoff = Double.parseDouble(manager.getSetting(SETTING.DE_FC_CUTOFF));
		deMinPctCutoff = Double.parseDouble(manager.getSetting(SETTING.DE_MIN_PCT_CUTOFF));
		netPvCutoff = Double.parseDouble(manager.getSetting(SETTING.NET_PV_CUTOFF));
		netLogFCCutoff = Double.parseDouble(manager.getSetting(SETTING.NET_FC_CUTOFF));
		positiveOnly = Boolean.parseBoolean(manager.getSetting(SETTING.POSITIVE_ONLY));
		heatMapCount = Integer.parseInt(manager.getSetting(SETTING.HEATMAP_COUNT));
		boolean viewOnly = Boolean.parseBoolean(manager.getSetting(SETTING.DONT_ANALYZE));
		if (viewOnly)
			doubleClickAction.setSelectedValue("View Data");
		else
			doubleClickAction.setSelectedValue("Create Network");
	}

	public void run(TaskMonitor monitor) {
		// Update all of our values
		manager.setSetting(SETTING.MAX_GENES, maxGenes);

		if (topGenes == 0)
			manager.setSetting(SETTING.TOP_GENES, "");
		else
			manager.setSetting(SETTING.TOP_GENES, topGenes);

		manager.setSetting(SETTING.DE_FC_CUTOFF, deLogFCCutoff);
		manager.setSetting(SETTING.DE_MIN_PCT_CUTOFF, deMinPctCutoff);
		manager.setSetting(SETTING.NET_PV_CUTOFF, netPvCutoff);
		manager.setSetting(SETTING.NET_FC_CUTOFF, netLogFCCutoff);
		if (doubleClickAction.getSelectedValue().equals("View Data"))
			manager.setSetting(SETTING.DONT_ANALYZE, "true");
		else
			manager.setSetting(SETTING.DONT_ANALYZE, "false");
		//
	}

	@ProvidesTitle
	public String getTitle() { return "scNetViz Settings"; }

}
