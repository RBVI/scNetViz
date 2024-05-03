package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;

import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTask;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.tasks.AbstractEmbeddingTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.CalculateDETask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ExportCSVTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteTSNETask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteUMAPTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteGraphTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.tSNETask;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;

public class CategoriesTab extends JPanel implements TaskObserver {
	final ScNVManager manager;
	final Experiment experiment;
	final ExperimentFrame expFrame;
	final CategoriesTab thisComponent;
	final List<Category> categoriesList;
	final List<String> categoriesNameList;
	final Map<Category, JTable> categoryTables;
	JTextField logFC;
	JTextField dDRThreshold;
	JComboBox<String> categories;
	JButton diffExpButton;
	public JButton cellPlotButton;
	final String accession;

	Category currentCategory = null;
	JScrollPane categoryPane;

	public CategoriesTab(final ScNVManager manager, final Experiment experiment, final ExperimentFrame expFrame) {
		this.manager = manager;
		this.experiment = experiment;
		this.categoriesList = experiment.getCategories();
		this.categoriesNameList = new ArrayList<>();
		this.categoryTables = new HashMap<>();
		for (Category c: categoriesList) {
			if (c != null)
				categoriesNameList.add(c.toString());
		}

		this.setLayout(new BorderLayout());
		thisComponent = this;	// Access to inner classes
		this.expFrame = expFrame;
		this.accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
		init();
	}

	public void selectAssays(List<String> assayList) {
		JTable catTable = categoryTables.get(currentCategory);

		// Clear the selection list
		catTable.clearSelection();
		catTable.setColumnSelectionAllowed(true);
		catTable.setRowSelectionAllowed(false);

		// Get the unsorted row labels
		List<String> colLabels = currentCategory.getMatrix().getColLabels(0);

		/*
		System.out.println("ColLabels = ");
		for (String lbl: colLabels) {
			System.out.println(lbl);
		}
		*/

		for (String assay: assayList) {
			// System.out.println("Selecting assay: "+assay);
			int index = colLabels.indexOf(assay);
			index = catTable.convertColumnIndexToView(index);
			catTable.setColumnSelectionInterval(index, index);
			catTable.scrollRectToVisible(new Rectangle(catTable.getCellRect(0, index, true)));
		}

		catTable.setColumnSelectionAllowed(false);
		catTable.setRowSelectionAllowed(true);
	}


	@Override
	public void allFinished(FinishStatus status) {
		diffExpButton.setEnabled(true);
	}

	@Override
	public void taskFinished(ObservableTask obsTask) {
		if (obsTask instanceof FileCategoryTask) {
			String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
			expFrame.addCategoriesContent(accession+": Categories Tab", new CategoriesTab(manager, experiment, expFrame));
		} else if (obsTask instanceof CalculateDETask) {
			DifferentialExpression diffExp = obsTask.getResults(DifferentialExpression.class);
			if (diffExp == null) {
				// System.out.println("diffExp = null!");
				SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						JOptionPane.showMessageDialog(expFrame, "Differential expression calculation failed", 
						                              "DE Failure", JOptionPane.ERROR_MESSAGE);
					}
				});
				diffExpButton.setEnabled(true);
				return;
			}
			DiffExpTab diffETab = new DiffExpTab(manager, experiment, expFrame, currentCategory, diffExp);
			expFrame.addDiffExpContent("Diff Exp", diffETab);
			diffExpButton.setEnabled(true);
		} else if (obsTask instanceof AbstractEmbeddingTask) {
			Map<String,double[]> embedding = ((AbstractEmbeddingTask)obsTask).getResults();
			if (embedding == null)
				return;
			experiment.setTSNE(embedding);
			showPlot();
			cellPlotButton.setEnabled(true);
			cellPlotButton.setText("View "+experiment.getPlotType());
			TPMTab tpmTab = expFrame.getTPMTab();
			if (tpmTab != null) {
				tpmTab.cellPlotButton.setEnabled(true);
				tpmTab.cellPlotButton.setText("View "+experiment.getPlotType());
			}
		}
		expFrame.toFront();
	}

	public void changeCategory(Category newCategory, int newRow) {
		// newCategory.setSelectedRow(newRow);
		JTable categoryTable = getCategoryTable(newCategory);
		categories.setSelectedItem(newCategory.toString());
		if (newRow >= 0)
			categoryTable.setRowSelectionInterval(newRow, newRow);
	}

	public void recalculateDE() {
		diffExpButton.doClick();
	}

	private void init() {
		JPanel buttonsPanelLeft = new JPanel();
		buttonsPanelLeft.setLayout(new BoxLayout(buttonsPanelLeft, BoxLayout.PAGE_AXIS));
		buttonsPanelLeft.setAlignmentX(Component.LEFT_ALIGNMENT);
		{
			//JLabel experimentLabel = new ExperimentLabel(experiment);
			//experimentLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
			//buttonsPanelLeft.add(experimentLabel);
			// buttonsPanelLeft.add(Box.createRigidArea(new Dimension(0, 10)));

			JLabel lbl = new JLabel("Available categories:");
			lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
			lbl.setAlignmentX(Component.LEFT_ALIGNMENT);
			buttonsPanelLeft.add(Box.createRigidArea(new Dimension(10,0)));
			buttonsPanelLeft.add(lbl);

			categories = new JComboBox<String>(categoriesNameList.toArray(new String[1]));
      categories.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					currentCategory = categoriesList.get(categories.getSelectedIndex());
					JTable table = getCategoryTable(currentCategory);
					SortableTableModel model = (SortableTableModel)table.getModel();
					categoryPane.setViewportView(table);
					// model.fireTableDataChanged();
					categoryPane.revalidate();
					categoryPane.repaint();
				}
			});
			categories.setFont(new Font("SansSerif", Font.PLAIN, 10));
      categories.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
				}
			});
			categories.setPreferredSize(new Dimension(200, 25));
			categories.setMaximumSize(new Dimension(200, 25));
			categories.setSize(new Dimension(200, 25));
			categories.setAlignmentX(Component.LEFT_ALIGNMENT);

			buttonsPanelLeft.add(Box.createRigidArea(new Dimension(10,0)));
			buttonsPanelLeft.add(categories);
		}

		JPanel centerPanel = new JPanel();
		centerPanel.setLayout(new BoxLayout(centerPanel, BoxLayout.PAGE_AXIS));
		{
			centerPanel.add(Box.createVerticalGlue());
			centerPanel.add(Box.createRigidArea(new Dimension(20,0)));
			{
				JLabel lbl = new JLabel("Cutoffs:");
				lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				lbl.setAlignmentX(Component.LEFT_ALIGNMENT);
				centerPanel.add(lbl);
			}

			JPanel settingsPanel = new JPanel();
			settingsPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
			settingsPanel.setLayout(new BoxLayout(settingsPanel, BoxLayout.LINE_AXIS));
			{
				settingsPanel.add(Box.createRigidArea(new Dimension(5,0)));
				JLabel log2FCLabel = new JLabel("Log2FC:");
				log2FCLabel.setFont(new Font("SansSerif", Font.BOLD, 10));
				log2FCLabel.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(log2FCLabel);
			}

			{
				logFC = new JTextField(manager.getSetting(SETTING.DE_FC_CUTOFF));
				logFC.setFont(new Font("SansSerif", Font.PLAIN, 10));
				logFC.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(logFC);
				settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));
			}

			{
				JLabel dDRThreshLabel = new JLabel("Min.pct:");
				dDRThreshLabel.setFont(new Font("SansSerif", Font.BOLD, 10));
				dDRThreshLabel.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(dDRThreshLabel);
			}

			{
				dDRThreshold = new JTextField(manager.getSetting(SETTING.DE_MIN_PCT_CUTOFF)+"%");
				dDRThreshold.setFont(new Font("SansSerif", Font.PLAIN, 10));
				dDRThreshold.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(dDRThreshold);
				settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));
			}

			{
				diffExpButton = new JButton("Calculate Diff Exp");
				diffExpButton.setFont(new Font("SansSerif", Font.BOLD, 10));
				diffExpButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						double log2FCCutoff = Double.parseDouble(logFC.getText());
						String dDRText = dDRThreshold.getText();
						double dDRCutoff = (double)Integer.parseInt(dDRText.replaceAll("%",""))/100.0;

						diffExpButton.setEnabled(false);
						TaskIterator ti = new TaskIterator(new CalculateDETask(manager, currentCategory, dDRCutoff, log2FCCutoff));
						try {
							manager.executeTasks(ti, thisComponent);
						} catch (Exception ex) {
							diffExpButton.setEnabled(true);
						}
					}
				});
				settingsPanel.add(diffExpButton);
			}

			centerPanel.add(settingsPanel);

		}
		
		JPanel buttonsPanelRight = new JPanel(new GridLayout(2,2));
		{
			JButton export = new JButton("Export CSV");
			export.setFont(new Font("SansSerif", Font.PLAIN, 10));
      export.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					ExportCSVTask task = new ExportCSVTask(manager, currentCategory.getMatrix());
					manager.executeTasks(new TaskIterator(task));
				}
			});
			
			cellPlotButton = new JButton("View Cell Plot");
			cellPlotButton.setFont(new Font("SansSerif", Font.PLAIN, 10));
			cellPlotButton.setEnabled(false);
      cellPlotButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					showPlot();
				}
			});

			buttonsPanelRight.add(ViewUtils.createPlotMenu(manager, experiment, thisComponent));
			buttonsPanelRight.add(cellPlotButton);
			buttonsPanelRight.add(ViewUtils.createCategoryMenu(manager, experiment));
			buttonsPanelRight.add(export);
		}

		JPanel topPanel = new JPanel(new BorderLayout());
		topPanel.add(buttonsPanelLeft, BorderLayout.WEST);
		topPanel.add(centerPanel, BorderLayout.CENTER);
		topPanel.add(buttonsPanelRight, BorderLayout.EAST);
		this.add(topPanel, BorderLayout.NORTH);

		JTable categoryTable = getCategoryTable(categoriesList.get(0));
		currentCategory = categoriesList.get(0);

		categoryPane = new JScrollPane(categoryTable);
		categoryPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		this.add(categoryPane, BorderLayout.CENTER);
		this.revalidate();
		this.repaint();
	}

	private JTable getCategoryTable(final Category category) {
		if (categoryTables.containsKey(category))
		 	return categoryTables.get(category);

		SortableTableModel tableModel = (SortableTableModel)category.getTableModel();
		// if (tableModel == null)
		// 	tableModel = new DefaultCategoryTableModel(category);

		JTable categoryTable = new SimpleTable(manager, tableModel);
		categoryTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent event) {
				int row = categoryTable.getSelectedRow();
				if (row >= 0) {
					// System.out.println("Setting selected row for "+category+" to "+row);
					category.setSelectedRow(categoryTable.convertRowIndexToModel(row));
				}
			}
		});
		categoryTables.put(category, categoryTable);
		return categoryTable;
	}

	private void showPlot() {
		String title = null;
		Category cat = currentCategory;
		int catRow = -1;
		if (cat != null) {
			title = experiment.getPlotType()+" Plot for "+accession+" Category "+cat.toString();
			catRow = cat.getSelectedRow();
			if (catRow < 0)
				catRow = cat.getDefaultRow();
			if (catRow >= 0) {
				List<String> rowLabels = cat.getMatrix().getRowLabels(0);
				title += " ("+rowLabels.get(catRow)+")";
			}
		}
		ViewUtils.showtSNE(manager, experiment, cat, catRow, -1, title);
	}

}
