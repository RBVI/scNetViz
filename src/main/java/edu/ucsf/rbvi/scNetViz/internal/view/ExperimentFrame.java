package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.cytoscape.util.swing.IconManager;

import org.cytoscape.work.TaskIterator;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.tasks.SettingsTask;

public class ExperimentFrame extends JFrame {
	final ScNVManager scManager;
	final String[] titles = {"TPM", "Categories", "DiffExp"};
	final JFrame jFrame;
	final Experiment experiment;
	final Font iconFont;
	JTabbedPane tabbedPane;
	JPanel headerPane;

	TPMTab tpmTab;
	CategoriesTab categoriesTab;
	DiffExpTab diffExpTab;

	public ExperimentFrame(final ScNVManager scManager, final Experiment experiment) {
		this.scManager = scManager;
		this.setLayout(new BorderLayout());
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		this.experiment = experiment;
		iconFont = scManager.getService(IconManager.class).getIconFont(17.0f);
		init();
		jFrame = this;
	}

	public void addTPMContent(String title, JPanel content) {
		titles[0] = title;
		jFrame.setTitle(title);
		JPanel panel = (JPanel)tabbedPane.getComponentAt(0);
		panel.removeAll();
		panel.add(content, BorderLayout.CENTER);
		tabbedPane.setEnabledAt(0, true);
		tpmTab = (TPMTab)content;
	}

	public void addCategoriesContent(String title, JPanel content) {
		titles[1] = title;
		JPanel panel = (JPanel)tabbedPane.getComponentAt(1);
		panel.removeAll();
		panel.add(content, BorderLayout.CENTER);
		tabbedPane.setEnabledAt(1, true);
		tabbedPane.setSelectedIndex(1);
		categoriesTab = (CategoriesTab)content;
	}

	public void addDiffExpContent(String title, JPanel content) {
		titles[2] = title;
		JPanel panel = (JPanel)tabbedPane.getComponentAt(2);
		panel.removeAll();
		if (content == null) {
			tabbedPane.setEnabledAt(2, false);
			diffExpTab = null;
			return;
		}
		panel.add(content, BorderLayout.CENTER);
		tabbedPane.setEnabledAt(2, true);
		tabbedPane.setSelectedIndex(2);
		diffExpTab = (DiffExpTab)content;
	}

	public TPMTab getTPMTab() { return tpmTab; }
	public CategoriesTab getCategoriesTab() { return categoriesTab; }
	public DiffExpTab getDiffExpTab() { return diffExpTab; }

	public void selectTab(int tab) {
		tabbedPane.setSelectedIndex(tab);
	}

	public void selectTab(String tab) {
		switch(tab) {
			case "tpmTable":
				selectTab(0);
				break;
			case "catTable":
				selectTab(1);
				break;
			case "DETable":
				selectTab(2);
				break;
		}
	}

	public void selectGenes(List<String> geneLabels) {
		tpmTab.selectGenes(geneLabels);
		diffExpTab.selectGenes(geneLabels);
	}

	public void selectAssays(List<String> assayLabels) {
		categoriesTab.selectAssays(assayLabels);
	}


	private void init() {
		headerPane = new JPanel();
		headerPane.setLayout(new BoxLayout(headerPane, BoxLayout.LINE_AXIS));
		headerPane.add(Box.createRigidArea(new Dimension(10,0)));
		{
			JButton helpButton = new HelpButton(scManager, null);
			headerPane.add(helpButton);
		}

		{
			JTextPane experimentLabel = new ExperimentLabel(experiment, headerPane.getBackground());
			headerPane.add(Box.createHorizontalGlue());
			headerPane.add(experimentLabel);
			headerPane.add(Box.createHorizontalGlue());
		}

		{
			JButton settings = new JButton(IconManager.ICON_COG);
			settings.setFont(iconFont);
			settings.setBorderPainted(false);
			settings.setContentAreaFilled(false);
			settings.setFocusPainted(false);
			settings.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					TaskIterator tasks = new TaskIterator(new SettingsTask(scManager));
					scManager.executeTasks(tasks);
				}
			});
			headerPane.add(settings);
			headerPane.add(Box.createRigidArea(new Dimension(10,0)));
		}

		this.add(headerPane, BorderLayout.NORTH);

		tabbedPane = new JTabbedPane();
		tabbedPane.setPreferredSize(new Dimension(1100, 500));
		tabbedPane.setFont(new Font("SansSerif", Font.BOLD, 10));
		this.add(tabbedPane, BorderLayout.CENTER);

		// OK, now add our three tabs with empty content
		tabbedPane.addTab("TPM", new EmptyJPanel());
		tabbedPane.addTab("Categories", new EmptyJPanel());
		tabbedPane.addTab("DiffExp", new EmptyJPanel());

		tabbedPane.setEnabledAt(0, false);
		tabbedPane.setEnabledAt(1, false);
		tabbedPane.setEnabledAt(2, false);

		// Set it up so the components provide the titles for the entire frame
		tabbedPane.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				int index = ((JTabbedPane)e.getSource()).getSelectedIndex();
				jFrame.setTitle(titles[index]);
			}
		});

		this.pack();
		this.setVisible(true);
		this.toFront();
	}

	class EmptyJPanel extends JPanel {
		EmptyJPanel() {
			super(new BorderLayout());
		}
	}
}
