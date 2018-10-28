package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ScrollPaneConstants;
import javax.swing.table.TableModel;

import java.util.HashMap;
import java.util.Map;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAExperiment;

public class TPMTab extends JPanel {
	final ScNVManager manager;
	final Experiment experiment;
	final TPMTab thisComponent;

	public TPMTab(final ScNVManager manager, final Experiment experiment) {
		this.manager = manager;
		this.experiment = experiment;

		this.setLayout(new BorderLayout());
		thisComponent = this;	// Access to inner classes
		init();
	}
	
	private void init() {

		JLabel experimentLabel = new ExperimentLabel(experiment);

		JPanel buttonsPanelRight = new JPanel(new GridLayout(1, 2));
		{
			JButton tsne = new JButton("View tSNE");
			tsne.setFont(new Font("SansSerif", Font.PLAIN, 10));
      tsne.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					if (experiment instanceof GXAExperiment) {
						String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
						String uri = "https://www.ebi.ac.uk/gxa/sc/experiments/"+accession+"/Results";
						Map<String, Object> args = new HashMap<>();
						args.put("newTab", "true");
						args.put("id", "GXA");
						args.put("url", uri);

						manager.executeCommand("cybrowser", "dialog", args);
					}
				}
			});
			buttonsPanelRight.add(tsne);
		}
		
		{
			JButton export = new JButton("Export CSV");
			export.setFont(new Font("SansSerif", Font.PLAIN, 10));
      export.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
				}
			});
			buttonsPanelRight.add(export);
		}
		
		JPanel topPanel = new JPanel(new BorderLayout());
		topPanel.add(experimentLabel, BorderLayout.WEST);
		topPanel.add(buttonsPanelRight, BorderLayout.EAST);
		this.add(topPanel, BorderLayout.NORTH);

		TableModel tableModel = experiment.getTableModel();
		if (tableModel == null)
			tableModel = new DefaultExperimentTableModel(experiment);

		JTable experimentTable = new ExperimentTable(manager, tableModel);

		JScrollPane scrollPane = new JScrollPane(experimentTable);
		scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		this.add(scrollPane, BorderLayout.CENTER);
		this.revalidate();
		this.repaint();
	}
}
