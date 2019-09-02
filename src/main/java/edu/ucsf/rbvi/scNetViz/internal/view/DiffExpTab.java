package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.BevelBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.tasks.CreateNetworkTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ExportCSVTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.HeatMapTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ViolinDiffExpTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ViolinGeneTask;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

public class DiffExpTab extends JPanel {
	final ScNVManager manager;
	final Experiment experiment;
	final ExperimentFrame expFrame;
	final List<Category> categories;
	final Category currentCategory;
	final DifferentialExpression diffExp;
	final DiffExpTab thisComponent;
	final Map<Category, List<String>> categoryLabelMap;
	JTable diffExpTable = null;

	JTextField topGenes;
	JTextField pValue;
	JTextField log2FC;
	JCheckBox positiveOnly;

	public DiffExpTab(final ScNVManager manager, final Experiment experiment, 
	                  final ExperimentFrame expFrame, final Category currentCategory,
										final DifferentialExpression diffExp) {
		this.manager = manager;
		this.experiment = experiment;
		this.setLayout(new BorderLayout());
		thisComponent = this;	// Access to inner classes
		this.expFrame = expFrame;
		this.categories = experiment.getCategories();
		this.currentCategory = currentCategory;
		this.diffExp = diffExp;
		categoryLabelMap = new HashMap<>();
		for (Category cat: categories) {
			categoryLabelMap.put(cat, cat.getMatrix().getRowLabels());
		}
		init();
	}

	public void selectGenes(List<String> geneList) {
		// Clear the selection list
		diffExpTable.clearSelection();
		diffExpTable.setRowSelectionAllowed(true);
		// Get the unsorted row labels
		List<String> rowLabels = experiment.getMatrix().getRowLabels();
		for (String gene: geneList) {
			int index = rowLabels.indexOf(gene);
			index = diffExpTable.convertRowIndexToView(index);
			diffExpTable.getSelectionModel().addSelectionInterval(index, index);
			diffExpTable.scrollRectToVisible(new Rectangle(diffExpTable.getCellRect(index, 0, true)));
		}
		String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
		ModelUtils.selectNodes(manager, accession, geneList);
	}

	private void init() {

		
		JPanel buttonsPanelLeft = new JPanel();
		buttonsPanelLeft.setLayout(new BoxLayout(buttonsPanelLeft, BoxLayout.PAGE_AXIS));
		{
			{
				JLabel lbl = new JLabel("Comparison:");
				lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				lbl.setAlignmentX(Component.LEFT_ALIGNMENT);
				buttonsPanelLeft.add(Box.createRigidArea(new Dimension(20,5)));
				buttonsPanelLeft.add(lbl);
			}

			{
				List<String> labels = new ArrayList<>();
				for (List<String> lbl: categoryLabelMap.values())
					labels.addAll(lbl);

				String selectedLabel = categoryLabelMap.get(currentCategory).get(diffExp.getCategoryRow());

				JComboBox<String> categoryBox = 
					new JComboBox<String>(labels.toArray(new String[1]));

				categoryBox.setSelectedItem(selectedLabel);

				// TODO: set our default value using currentCategory and the selectedRow
 	     	categoryBox.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						String newLabel = (String)categoryBox.getSelectedItem();

						for (Category cat: categoryLabelMap.keySet()) {
							for (int i = 0; i < categoryLabelMap.get(cat).size(); i++) {
								if (newLabel.equals(categoryLabelMap.get(cat).get(i))) {
									changeCategory(cat, i);
									return;
								}
							}
						}
					}

					void changeCategory(Category cat, int selectedRow) {
						CategoriesTab catTab = expFrame.getCategoriesTab();
						catTab.changeCategory(cat, selectedRow);
						// Recalculate
						catTab.recalculateDE();
						// Refresh
					}
				});
			 	Dimension size = new Dimension(250,25);
				categoryBox.setPreferredSize(size);
				categoryBox.setMaximumSize(size);
				categoryBox.setSize(size);
				categoryBox.setFont(new Font("SansSerif", Font.PLAIN, 10));
				categoryBox.setAlignmentX(Component.LEFT_ALIGNMENT);
				buttonsPanelLeft.add(categoryBox);
			}

			{
				JComboBox comparison = new JComboBox(getComparisons());
				comparison.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
					}
				});
			 	Dimension size = new Dimension(250,25);
				comparison.setPreferredSize(size);
				comparison.setMaximumSize(size);
				comparison.setSize(size);
				comparison.setFont(new Font("SansSerif", Font.PLAIN, 10));
				comparison.setAlignmentX(Component.LEFT_ALIGNMENT);
				buttonsPanelLeft.add(comparison);
			}
		}

		JPanel centerPanel = new JPanel();
		centerPanel.setLayout(new BoxLayout(centerPanel, BoxLayout.PAGE_AXIS));
		{
			centerPanel.add(Box.createVerticalGlue());
			centerPanel.add(Box.createRigidArea(new Dimension(10,0)));
			{
				JLabel lbl = new JLabel("Network analysis:");
				lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				lbl.setAlignmentX(Component.LEFT_ALIGNMENT);
				centerPanel.add(lbl);
			}

			JPanel settingsPanel = new JPanel();
			settingsPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
			settingsPanel.setLayout(new BoxLayout(settingsPanel, BoxLayout.LINE_AXIS));

			JPanel settingsPanel1 = new JPanel();
			settingsPanel1.setAlignmentX(Component.LEFT_ALIGNMENT);
			settingsPanel1.setLayout(new BoxLayout(settingsPanel1, BoxLayout.LINE_AXIS));
			settingsPanel1.setBorder(BorderFactory.createEtchedBorder());
			{
				settingsPanel1.add(Box.createRigidArea(new Dimension(5,0)));
				JLabel pValueLbl = new JLabel("pValue:");
				pValueLbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				pValueLbl.setMaximumSize(new Dimension(50,35));
				settingsPanel1.add(pValueLbl);
			}

			{
				pValue = new JTextField(manager.getSetting(SETTING.NET_PV_CUTOFF));
				pValue.setFont(new Font("SansSerif", Font.PLAIN, 10));
				pValue.setMaximumSize(new Dimension(50,35));
				settingsPanel1.add(pValue);
				settingsPanel1.add(Box.createRigidArea(new Dimension(15,0)));
			}

			{
				JLabel log2FCLbl = new JLabel("Log2FC:");
				log2FCLbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				log2FCLbl.setMaximumSize(new Dimension(50,35));
				settingsPanel1.add(log2FCLbl);
			}

			{
				log2FC = new JTextField(manager.getSetting(SETTING.NET_FC_CUTOFF));
				log2FC.setFont(new Font("SansSerif", Font.PLAIN, 10));
				log2FC.setMaximumSize(new Dimension(50,35));
				settingsPanel1.add(log2FC);
				settingsPanel1.add(Box.createRigidArea(new Dimension(5,0)));
			}

			settingsPanel.add(settingsPanel1);
			settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));

			{
				JLabel orLbl = new JLabel(" OR ");
				orLbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				orLbl.setMaximumSize(new Dimension(35,35));
				settingsPanel.add(orLbl);
				settingsPanel.add(Box.createRigidArea(new Dimension(5,0)));
			}

			JPanel settingsPanel2 = new JPanel();
			settingsPanel2.setAlignmentX(Component.LEFT_ALIGNMENT);
			settingsPanel2.setLayout(new BoxLayout(settingsPanel2, BoxLayout.LINE_AXIS));
			settingsPanel2.setBorder(BorderFactory.createEtchedBorder());
			settingsPanel2.add(Box.createRigidArea(new Dimension(5,0)));
			{
				JLabel topGenesLbl = new JLabel("Top n genes");
				topGenesLbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				topGenesLbl.setMaximumSize(new Dimension(70,35));
				settingsPanel2.add(topGenesLbl);
			}

			{
				topGenes = new JTextField(manager.getSetting(SETTING.TOP_GENES));
				topGenes.setFont(new Font("SansSerif", Font.PLAIN, 10));
				topGenes.setMaximumSize(new Dimension(50,35));
				settingsPanel2.add(topGenes);
			}

			settingsPanel.add(settingsPanel2);
			settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));

			{
				positiveOnly = new JCheckBox("Positive only", 
				                             Boolean.parseBoolean(manager.getSetting(SETTING.POSITIVE_ONLY)));
				positiveOnly.setFont(new Font("SansSerif", Font.BOLD, 10));
				settingsPanel.add(positiveOnly);
				settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));
			}

			{
				JButton createNetwork = new JButton("Create Network");
				createNetwork.setFont(new Font("SansSerif", Font.BOLD, 10));
				createNetwork.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						int topNGenes = -1;
						double log2FCCutoff = 1.0;
						double pv = .05;

						if (topGenes.getText().length() > 0)
							topNGenes = Integer.parseInt(topGenes.getText());
						if (log2FC.getText().length() > 0)
							log2FCCutoff = Double.parseDouble(log2FC.getText());
						if (pValue.getText().length() > 0)
							pv = Double.parseDouble(pValue.getText());

						int maxGenes = Integer.parseInt(manager.getSetting(SETTING.MAX_GENES));
						boolean posOnly = Boolean.parseBoolean(manager.getSetting(SETTING.POSITIVE_ONLY));

						createNetwork.setEnabled(false);

						// TODO: Add max genes somewhere
						// TODO: Get the result of the network creation to support selection
						posOnly = positiveOnly.isSelected();
						CreateNetworkTask task = new CreateNetworkTask(manager, diffExp, pv, log2FCCutoff, topNGenes,
						                                               posOnly, maxGenes);
						manager.executeTasks(new TaskIterator(task));
						createNetwork.setEnabled(true);
					}
				});
				settingsPanel.add(createNetwork);
			}

			centerPanel.add(settingsPanel);

		}

		JPanel buttonsPanelRight = new JPanel(/*new GridLayout(3, 1)*/);
		buttonsPanelRight.setLayout(new BoxLayout(buttonsPanelRight, BoxLayout.PAGE_AXIS));
		{
			/*
			JButton viewHeatMap = new JButton("View Heatmap");
			viewHeatMap.setFont(new Font("SansSerif", Font.PLAIN, 10));
      viewHeatMap.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					Map<String, double[]> dataMap = new HashMap<>();
					List<String> columns = new ArrayList<>();
					boolean posOnly = positiveOnly.isSelected();
					for (Object cat: diffExp.getLogGERMap().keySet()) {
						double[] logGER = diffExp.getLogGER(cat, posOnly);
						if (logGER != null) {
							dataMap.put(currentCategory.mkLabel(cat), logGER);
							columns.add(currentCategory.mkLabel(cat));
						}
					}
					// Use a separate task for this since we've got some options...
					HeatMapTask task = new HeatMapTask(manager, currentCategory, diffExp.getRowLabels(), 
					                                   dataMap, columns, posOnly, -1, null);
					manager.executeTasks(new TaskIterator(task));
				}
			});
			*/
			buttonsPanelRight.add(new JLabel(" "));

			{
				Map<String, Task> map = new LinkedHashMap<>();
				map.put("Heatmap", new HeatMapTaskWrapper());
				map.put("Violin (diff exp)", new ViolinDiffExpTaskWrapper());
				map.put("Violin (gene)", new ViolinGeneTaskWrapper());
				buttonsPanelRight.add(new PullDownMenu(manager, "View Plots", map, null));
			}

			{
				JButton export = new JButton("Export CSV");
				export.setFont(new Font("SansSerif", Font.PLAIN, 10));
 	  	   export.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						ExportCSVTask task = new ExportCSVTask(manager, diffExp);
						manager.executeTasks(new TaskIterator(task));
					}
				});
				buttonsPanelRight.add(export);
			}

			/*
			JButton viewViolin = new JButton("View Violin Plot");
			viewViolin.setFont(new Font("SansSerif", Font.PLAIN, 10));
      viewViolin.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					Map<String, double[]> dataMap = new HashMap<>();
					List<String> columns = new ArrayList<>();
					for (Object cat: diffExp.getLogGERMap().keySet()) {
						double[] logGER = diffExp.getLogGER(cat, positiveOnly.isSelected());
						if (logGER != null) {
							dataMap.put(currentCategory.mkLabel(cat), logGER);
							columns.add(currentCategory.mkLabel(cat));
						}
					}
					String[] dataAndNames = CyPlotUtils.mapToDataAndNames(dataMap, diffExp.getRowLabels(), columns);
					String accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
					String title = experiment.getSource().toString()+" "+ accession+ " Differential Expression";
					CyPlotUtils.createViolinPlot(manager, dataAndNames[0], dataAndNames[1], 
					                             CyPlotUtils.listToCSV(columns), title, "", "Log(FC)", accession);
				}
			});

			buttonsPanelRight.add(viewHeatMap);
			buttonsPanelRight.add(viewViolin);
			buttonsPanelRight.add(export);
			*/
		}

		JPanel topPanel = new JPanel(new BorderLayout());
		topPanel.add(buttonsPanelLeft, BorderLayout.WEST);
		topPanel.add(centerPanel, BorderLayout.CENTER);
		topPanel.add(buttonsPanelRight, BorderLayout.EAST);
		this.add(topPanel, BorderLayout.NORTH);

		SortableTableModel tableModel = diffExp.getTableModel();

		diffExpTable = new SimpleTable(manager, tableModel);
		diffExpTable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		diffExpTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent event) {
				int[] rows = diffExpTable.getSelectedRows();
				if (rows.length == 0) return;
				List<String> geneList = new ArrayList<>();
				for (int row: rows) {
					geneList.add(diffExpTable.getValueAt(row, 0).toString());
				}
				String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
				ModelUtils.selectNodes(manager, accession, geneList);
			}
		});

		JScrollPane diffExpPane = new JScrollPane(diffExpTable);
		diffExpPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		this.add(diffExpPane, BorderLayout.CENTER);

		this.revalidate();
		this.repaint();
	}

	// TODO: This should be more dynamic
	String[] getComparisons() {
		String[] comparisonTypes = {"Each vs. Others"};
		return comparisonTypes;
	}

	class HeatMapTaskWrapper extends AbstractTask {
		public void run(TaskMonitor monitor) {
			boolean posOnly = positiveOnly.isSelected();
			Task heatTask = new HeatMapTask(manager, currentCategory, diffExp, posOnly, -1, null);
			insertTasksAfterCurrentTask(heatTask);
		}
	}

	class ViolinDiffExpTaskWrapper extends AbstractTask {
		public void run(TaskMonitor monitor) {
			boolean posOnly = positiveOnly.isSelected();
			Task vdeTask = new ViolinDiffExpTask(manager, currentCategory, diffExp, posOnly);
			insertTasksAfterCurrentTask(vdeTask);
		}
	}

	class ViolinGeneTaskWrapper extends AbstractTask {
		public void run(TaskMonitor monitor) {
			int selectedCategory = diffExp.getCategoryRow();
			int[] rows = diffExpTable.getSelectedRows();
			if (rows == null || rows.length == 0) {
				SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						JOptionPane.showMessageDialog(expFrame, "Gene row must be selected!", 
						                              "No gene", JOptionPane.ERROR_MESSAGE);
					}

				});
				return;
			}
			Task vgeneTask = new ViolinGeneTask(manager, currentCategory, selectedCategory, rows[0]);
			insertTasksAfterCurrentTask(vgeneTask);
		}
	}
}
