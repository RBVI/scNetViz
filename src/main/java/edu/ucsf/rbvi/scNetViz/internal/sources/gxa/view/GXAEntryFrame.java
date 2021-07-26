package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.view;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

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
import javax.swing.JTextField;
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

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.model.Species;
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
			// Right panel == settings & search
			// Middle panel == two buttons
			JPanel rightPanel = new JPanel();
			rightPanel.setLayout(new BoxLayout(rightPanel, BoxLayout.LINE_AXIS));
			{
				JLabel label = new JLabel(IconManager.ICON_SEARCH);
				label.setFont(iconFont);
				JTextField search = new JTextField(16);
				{
					search.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent e) {
							String searchText = search.getText();
							GXAEntryTableModel model = (GXAEntryTableModel)gxaEntryTable.getModel();
							gxaEntryTable.clearSelection();
							List<Integer> matches = model.search(searchText);
							for (Integer row: matches) {
								int modelRow = gxaEntryTable.getRowSorter().convertRowIndexToView(row);
								gxaEntryTable.getSelectionModel().addSelectionInterval(modelRow, modelRow);
                gxaEntryTable.scrollRectToVisible(gxaEntryTable.getCellRect(modelRow, 0, true));
							}
							search.setText("");
						}
					});
					search.addFocusListener(new FocusListener() {
						public void focusGained(FocusEvent e) {
							search.setText("");
						}
						public void focusLost(FocusEvent e) {}
					});
				}
				search.setText("Search");
				search.setMaximumSize(new Dimension(200,30));
				rightPanel.add(label);
				rightPanel.add(search);
				rightPanel.add(Box.createRigidArea(new Dimension(10,0)));

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
						args.put("url", "http://www.rbvi.ucsf.edu/cytoscape/scNetViz/index.shtml");
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
							loadExperiment(selectedAcc, true);
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

			GXAEntryTableModel tableModel = new GXAEntryTableModel(gxaSource.getMetadata());
			gxaEntryTable = new GXAEntryTable(scNVManager, gxaSource, tableModel);
			scrollPane = new JScrollPane(gxaEntryTable);
			scrollPane.setPreferredSize(new Dimension(900, 500));
			this.add(scrollPane, BorderLayout.CENTER);
			this.pack();
			this.setVisible(true);
			this.toFront();
		}
	}

	public void enableButtons(String acc, boolean enable) {
		this.selectedAcc = acc;
    Metadata mdata = gxaSource.getMetadata(acc);
    if (acc == null) return;

    // See if we know the species for this accession
    if (Species.getSpecies(mdata.get(Metadata.SPECIES).toString()) != null)
		  createButton.setEnabled(enable);
    else
		  createButton.setEnabled(false);
		viewButton.setEnabled(enable);
	}

	public void loadExperiment(String acc) {
		boolean dontAnalyze = Boolean.parseBoolean(scNVManager.getSetting(SETTING.DONT_ANALYZE));
		loadExperiment(acc, dontAnalyze);
	}

	private void loadExperiment(String acc, boolean dontAnalyze) {
		// Have we already loaded this experiment
		Experiment exp = scNVManager.getExperiment(acc);
		if (exp == null) {
			// No, load it
			TaskIterator tasks = new TaskIterator(new LoadExperimentTask(acc, dontAnalyze));
			scNVManager.executeTasks(tasks);
		} else {
			TaskIterator tasks = new TaskIterator(new ShowExperimentTableTask(scNVManager, exp));
			scNVManager.executeTasks(tasks);
		}
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
			GXAExperiment experiment = (GXAExperiment)gxaSource.getExperiment(accession, taskMonitor, true);
			if (experiment != null)
				scNVManager.addExperiment(accession,experiment);
			gxaSource.showEntriesTable(false);

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
