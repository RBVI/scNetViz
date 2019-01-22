package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.util.Collections;
import java.util.List;
import java.util.Properties;
import javax.swing.SwingUtilities;

import org.cytoscape.application.events.SetCurrentNetworkListener;
import org.cytoscape.application.swing.CySwingApplication;
import org.cytoscape.application.swing.CytoPanel;
import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelComponent2;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.application.swing.CytoPanelState;
import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;
import edu.ucsf.rbvi.scNetViz.internal.view.CategoriesTab;
import edu.ucsf.rbvi.scNetViz.internal.view.DiffExpTab;
import edu.ucsf.rbvi.scNetViz.internal.view.ScNVCytoPanel;
import edu.ucsf.rbvi.scNetViz.internal.view.TPMTab;

public class ShowResultsPanelTask extends AbstractTask {
	final ScNVManager manager;
	Experiment experiment = null;

	@Tunable (description="Experiment to show")
	public ListSingleSelection<Experiment> accession = null;

	public ShowResultsPanelTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		accession = new ListSingleSelection<>(manager.getExperiments());
	}

	public ShowResultsPanelTask(final ScNVManager manager, Experiment experiment) {
		super();
		this.manager = manager;
		this.experiment = experiment;
		accession = new ListSingleSelection<Experiment>(Collections.singletonList(experiment));
	}

	public void run(TaskMonitor monitor) {
		if (accession != null)
			experiment = accession.getSelectedValue();

		ScNVCytoPanel resultsPanel = new ScNVCytoPanel(manager, experiment);
		manager.registerService(resultsPanel, CytoPanelComponent.class, new Properties());
		manager.registerService(resultsPanel, SetCurrentNetworkListener.class, new Properties());

		CySwingApplication swingApp = manager.getService(CySwingApplication.class);
		CytoPanel panel = swingApp.getCytoPanel(CytoPanelName.EAST);
		panel.setState(CytoPanelState.DOCK);
	}
}
