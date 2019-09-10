package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cytoscape.event.CyEventHelper;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.TaskObserver;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.json.JSONResult;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

// Tunable to choose experiment?

public class CreateNetworkCommandTask extends AbstractTask {
	final ScNVManager manager;
	final CyEventHelper cyEventHelper;
	CyNetwork unionNetwork = null;
	VisualStyle baseStyle = null;

	DifferentialExpression diffExp = null;

	@Tunable (description="Experiment accession")
	public ListSingleSelection<Experiment> accession;

	@Tunable (description="FDR cutoff")
	public double fdrCutoff;

	@Tunable (description="Log2FC cutoff")
	public double log2FCCutoff;

	@Tunable (description="The maximum number of genes")
	public int maxGenes;

	@Tunable (description="Only consider positive fold changes")
	public boolean positiveOnly;

	public CreateNetworkCommandTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		cyEventHelper = manager.getService(CyEventHelper.class);
		accession = new ListSingleSelection<>(manager.getExperiments());
		fdrCutoff = Double.parseDouble(manager.getSetting(SETTING.NET_PV_CUTOFF));
		log2FCCutoff = Double.parseDouble(manager.getSetting(SETTING.NET_FC_CUTOFF));
		maxGenes = Integer.parseInt(manager.getSetting(SETTING.MAX_GENES));
		positiveOnly = Boolean.parseBoolean(manager.getSetting(SETTING.POSITIVE_ONLY));
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("Creating Networks");
		if (accession != null) {
			diffExp = accession.getSelectedValue().getDiffExp();
			if (diffExp == null) {
				monitor.showMessage(TaskMonitor.Level.ERROR, 
				                    "No differential expression has been calculated for "+accession.getSelectedValue());
				return;
			}
		}

		insertTasksAfterCurrentTask(new CreateNetworkTask(manager, diffExp, fdrCutoff, log2FCCutoff, maxGenes, positiveOnly, maxGenes));
	}
}
