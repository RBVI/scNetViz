package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class ExperimentFrame extends JFrame {
	final ScNVManager scManager;
	final String[] titles = {"TPM", "Categories", "DiffExp"};
	final JFrame jFrame;
	JTabbedPane tabbedPane;

	TPMTab tpmTab;
	CategoriesTab categoriesTab;
	DiffExpTab diffExpTab;

	public ExperimentFrame(final ScNVManager scManager) {
		this.scManager = scManager;
		this.setLayout(new BorderLayout());
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
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
		panel.add(content, BorderLayout.CENTER);
		tabbedPane.setEnabledAt(2, true);
		tabbedPane.setSelectedIndex(2);
		diffExpTab = (DiffExpTab)content;
	}

	public TPMTab getTPMTab() { return tpmTab; }
	public CategoriesTab getCategoriesTab() { return categoriesTab; }
	public DiffExpTab getDiffExpTab() { return diffExpTab; }

	private void init() {
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
	}

	class EmptyJPanel extends JPanel {
		EmptyJPanel() {
			super(new BorderLayout());
		}
	}
}
