package edu.ucsf.rbvi.scNetViz.internal.tasks;

import javax.swing.SwingUtilities;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;
import edu.ucsf.rbvi.scNetViz.internal.view.TPMTab;
import edu.ucsf.rbvi.scNetViz.internal.view.CategoriesTab;

// Tunable to choose experiment?

public class ShowExperimentTableTask extends AbstractTask {
	final ScNVManager manager;
	Experiment experiment = null;

	public ShowExperimentTableTask(final ScNVManager manager) {
		super();
		this.manager = manager;
	}

	public ShowExperimentTableTask(final ScNVManager manager, Experiment experiment) {
		super();
		this.manager = manager;
		this.experiment = experiment;
	}

	public void run(TaskMonitor monitor) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				// Create the Experiment Frame
				ExperimentFrame frame = new ExperimentFrame(manager);
				// Add our TPM tab
				String accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
				System.out.println("Accession = "+accession);
				frame.addTPMContent(accession+": TPM Tab", new TPMTab(manager, experiment));
				// Add our Categories tab
				frame.addCategoriesContent(accession+": Categories Tab", new CategoriesTab(manager, experiment));
			}
		});
	}
}
