package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;

import javax.swing.JPanel;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class CategoriesTab extends JPanel {
	final ScNVManager manager;
	final Experiment experiment;
	final CategoriesTab thisComponent;

	public CategoriesTab(final ScNVManager manager, final Experiment experiment) {
		this.manager = manager;
		this.experiment = experiment;

		this.setLayout(new BorderLayout());
		thisComponent = this;	// Access to inner classes
	}
}
