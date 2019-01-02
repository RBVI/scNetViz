package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.view;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.io.File;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import org.cytoscape.application.swing.CytoPanelComponent2;
import org.cytoscape.application.swing.CytoPanelName;

import org.cytoscape.util.swing.IconManager;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAExperiment;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ProcessAllTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.SettingsTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ShowExperimentTableTask;

public class GXAEntryFrame extends JFrame {
	final ScNVManager scNVManager;
	JScrollPane scrollPane = null;
	GXAEntryTable gxaEntryTable;
	GXASource gxaSource;
	// JCheckBox loadOnly;
	Icon ico = null;
	JButton viewButton = null;
	JButton createButton = null;
	String selectedAcc = null;
	final Font iconFont;

	public GXAEntryFrame (final ScNVManager scNVManager, final GXASource gxaSource) {
		super();
		this.scNVManager = scNVManager;
		this.gxaSource = gxaSource;
		this.setLayout(new BorderLayout());
		this.setTitle("Single Cell Expression Atlas Browser");
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		iconFont = scNVManager.getService(IconManager.class).getIconFont(17.0f);
		init();
	}

	public void init() {
		if (scrollPane == null) {
			if (gxaSource.getMetadata().size() == 0) {
				gxaSource.loadGXAEntries(null);
			}

			// Left panel == "Select an experiment"
			// Right panel == settings
			// Middle panel == two buttons
			JPanel rightPanel = new JPanel();
			rightPanel.setLayout(new BoxLayout(rightPanel, BoxLayout.LINE_AXIS));
			{
				JButton settings = new JButton(IconManager.ICON_COG);
				settings.setFont(iconFont);
				settings.setBorderPainted(false);
				settings.setContentAreaFilled(false);
				settings.setFocusPainted(false);
				settings.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						TaskIterator tasks = new TaskIterator(new SettingsTask(scNVManager));
						scNVManager.executeTasks(tasks);
					}
				});
				rightPanel.add(settings);
				rightPanel.add(Box.createRigidArea(new Dimension(10,0)));
			}

			JPanel leftPanel = new JPanel();
			leftPanel.setLayout(new BoxLayout(leftPanel, BoxLayout.LINE_AXIS));
			{
				leftPanel.add(Box.createRigidArea(new Dimension(10,0)));

				JButton helpButton = new JButton("Help");
				helpButton.setFont(new Font("SansSerif", Font.PLAIN, 10));
				helpButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						Map<String, Object> args = new HashMap<>();
						args.put("id","scNetViz");
						args.put("title", "scNetViz Help");
						args.put("url", "http://preview.rbvi.ucsf.edu/cytoscape/scNetViz/index.shtml");
						scNVManager.executeCommand("cybrowser", "dialog", args, false);
					}
				});
				leftPanel.add(helpButton);
			}

			JPanel middlePanel = new JPanel();
			middlePanel.setLayout(new BoxLayout(middlePanel, BoxLayout.PAGE_AXIS));
			{
				middlePanel.add(Box.createHorizontalGlue());

				JLabel labl = new JLabel("Select an experiment to analyze");
				labl.setFont(new Font("SansSerif", Font.ITALIC, 12));
				labl.setAlignmentX(Component.CENTER_ALIGNMENT);
				middlePanel.add(labl);

				JPanel middleButtonPanel = new JPanel();
				middleButtonPanel.setLayout(new BoxLayout(middleButtonPanel, BoxLayout.LINE_AXIS));
				{
					middleButtonPanel.add(Box.createHorizontalGlue());
					viewButton = new JButton("View Data");
					viewButton.setFont(new Font("SansSerif", Font.BOLD, 10));
					viewButton.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent e) {
							TaskIterator tasks = new TaskIterator(new LoadExperimentTask(selectedAcc, true));
							scNVManager.executeTasks(tasks);
						}
					});
					viewButton.setEnabled(false);
					middleButtonPanel.add(viewButton);
	
					middleButtonPanel.add(Box.createRigidArea(new Dimension(10,0)));
	
					createButton = new JButton("Create Networks");
					createButton.setFont(new Font("SansSerif", Font.BOLD, 10));
					createButton.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent e) {
							TaskIterator tasks = new TaskIterator(new LoadExperimentTask(selectedAcc, false));
							scNVManager.executeTasks(tasks);
						}
					});
					createButton.setEnabled(false);
					middleButtonPanel.add(createButton);
					middleButtonPanel.add(Box.createHorizontalGlue());
				}

				middlePanel.add(middleButtonPanel);
				middlePanel.add(Box.createHorizontalGlue());
			}

			JPanel topPanel = new JPanel(new BorderLayout());
			topPanel.add(leftPanel, BorderLayout.WEST);
			topPanel.add(middlePanel, BorderLayout.CENTER);
			topPanel.add(rightPanel, BorderLayout.EAST);

			this.add(topPanel, BorderLayout.NORTH);


			/*
			boolean dontAnalyze = Boolean.parseBoolean(scNVManager.getSetting(SETTING.DONT_ANALYZE));

			JPanel bottomPanel = new JPanel(new BorderLayout());
			bottomPanel.add(new JLabel("<html><i>Double click on experiment row to launch network analysis</i></html>", SwingConstants.LEFT), 
			             BorderLayout.WEST);
			loadOnly = new JCheckBox("<html><span style=\"font-size: 80%\">Load experiment only (don't analyze)</span></html>");
			loadOnly.setHorizontalTextPosition(SwingConstants.LEFT);
			loadOnly.setSelected(dontAnalyze);
			bottomPanel.add(loadOnly, BorderLayout.EAST);

			this.add(bottomPanel, BorderLayout.SOUTH);
			*/

			GXAEntryTableModel tableModel = new GXAEntryTableModel(gxaSource.getMetadata());
			gxaEntryTable = new GXAEntryTable(scNVManager, gxaSource, tableModel);
			scrollPane = new JScrollPane(gxaEntryTable);
			scrollPane.setPreferredSize(new Dimension(900, 500));
			this.add(scrollPane, BorderLayout.CENTER);
			this.pack();
			this.setVisible(true);
		}
	}

	public void enableButtons(String acc, boolean enable) {
		this.selectedAcc = acc;
		viewButton.setEnabled(enable);
		createButton.setEnabled(enable);
	}

	public void loadExperiment(String acc) {
		boolean dontAnalyze = Boolean.parseBoolean(scNVManager.getSetting(SETTING.DONT_ANALYZE));
		TaskIterator tasks = new TaskIterator(new LoadExperimentTask(acc, dontAnalyze));
		scNVManager.executeTasks(tasks);
	}

	// public boolean isLoadOnly() {
// 		return loadOnly.isSelected();
	// }

	@Override
	public String getTitle() {
		return "Single Cell Expression Atlas Browser";
	}

	public class LoadExperimentTask extends AbstractTask {
		String accession;
		boolean isLoadOnly;

		LoadExperimentTask(String acc, boolean isLoadOnly) {
			this.accession = acc;
			this.isLoadOnly = isLoadOnly;
		}

		@Override
		public void run(TaskMonitor taskMonitor) {
			taskMonitor.setTitle("Loading Single Cell Expression Atlas experiment "+accession);
			GXAExperiment experiment = (GXAExperiment)gxaSource.getExperiment(accession, taskMonitor);
			if (experiment != null)
				scNVManager.addExperiment(accession,experiment);
			gxaSource.showEntriesTable(false);
			experiment.fetchClusters(taskMonitor);
			experiment.fetchDesign(taskMonitor);

			// Now, depending on the setting of our checkbox, either show
			// the experiments table or process the data
			if (isLoadOnly) {
				// Create the Experiment table and show it
				insertTasksAfterCurrentTask(new ShowExperimentTableTask(scNVManager, experiment));
			} else {
				// Create the autoprocess task and append it
				insertTasksAfterCurrentTask(new ProcessAllTask(scNVManager, experiment));
			}
		}

		@ProvidesTitle
		public String getTitle() {return "Load Single Cell Expression Atlas Experiment "+accession;}
	}

}
