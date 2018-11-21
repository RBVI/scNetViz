package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
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
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ScrollPaneConstants;
import javax.swing.table.TableModel;

import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTask;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTaskFactory;

public class CategoriesTab extends JPanel implements TaskObserver {
	final ScNVManager manager;
	final Experiment experiment;
	final ExperimentFrame expFrame;
	final CategoriesTab thisComponent;
	final List<Category> categoriesList;
	final List<String> categoriesNameList;
	final Map<Category, JTable> categoryTables;

	Category currentCategory = null;
	JScrollPane categoryPane;

	public CategoriesTab(final ScNVManager manager, final Experiment experiment, final ExperimentFrame expFrame) {
		this.manager = manager;
		this.experiment = experiment;
		this.categoriesList = experiment.getCategories();
		this.categoriesNameList = new ArrayList<>();
		this.categoryTables = new HashMap<>();
		for (Category c: categoriesList) {
			categoriesNameList.add(c.toString());
		}

		this.setLayout(new BorderLayout());
		thisComponent = this;	// Access to inner classes
		this.expFrame = expFrame;
		init();
	}

	@Override
	public void allFinished(FinishStatus status) {
	}

	@Override
	public void taskFinished(ObservableTask obsTask) {
		if (obsTask instanceof FileCategoryTask) {
			String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
			expFrame.addCategoriesContent(accession+": Categories Tab", new CategoriesTab(manager, experiment, expFrame));
		}
	}
	

	private void init() {

		JPanel buttonsPanelRight = new JPanel(new GridLayout(4,1));
		// buttonsPanelRight.setLayout(new BoxLayout(buttonsPanelRight, BoxLayout.PAGE_AXIS));
		
		{
			JButton importCategory = new JButton("Add Category");
			importCategory.setFont(new Font("SansSerif", Font.PLAIN, 10));
      importCategory.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					// We need to use the File importer for this
					TaskFactory importCategory = 
									new FileCategoryTaskFactory(manager, (FileSource)manager.getSource("file"), experiment);
					manager.executeTasks(importCategory, thisComponent);
				}
			});
			// buttonsPanelRight.add(importCategory);

			// JLabel lbl = new JLabel("");
			JButton export = new JButton("Export CSV");
			export.setFont(new Font("SansSerif", Font.PLAIN, 10));
      export.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
				}
			});
			JButton diffExp = new JButton("Calculate Diff Exp");
			diffExp.setFont(new Font("SansSerif", Font.PLAIN, 10));
      diffExp.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					DiffExpTab diffETab = new DiffExpTab(manager, experiment, expFrame, currentCategory);
					expFrame.addDiffExpContent("Diff Exp", diffETab);
				}
			});

			// buttonsPanelRight.add(lbl);
			// buttonsPanelRight.add(Box.createVerticalGlue());
			buttonsPanelRight.add(new JLabel(""));
			buttonsPanelRight.add(importCategory);
			buttonsPanelRight.add(export);
			buttonsPanelRight.add(diffExp);
		}
		
		JPanel buttonsPanelLeft = new JPanel();
		buttonsPanelLeft.setLayout(new BoxLayout(buttonsPanelLeft, BoxLayout.PAGE_AXIS));
		buttonsPanelLeft.setAlignmentX(Component.LEFT_ALIGNMENT);
		{
			JLabel experimentLabel = new ExperimentLabel(experiment);
			experimentLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
			buttonsPanelLeft.add(experimentLabel);
			buttonsPanelLeft.add(Box.createRigidArea(new Dimension(0, 10)));

			JLabel lbl = new JLabel("Available categories:");
			lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
			lbl.setAlignmentX(Component.LEFT_ALIGNMENT);
			buttonsPanelLeft.add(lbl);

			JComboBox categories = new JComboBox(categoriesNameList.toArray());
      categories.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					currentCategory = categoriesList.get(categories.getSelectedIndex());
					JTable table = getCategoryTable(currentCategory);
					SortableTableModel model = (SortableTableModel)table.getModel();
					categoryPane.setViewportView(table);
					model.fireTableDataChanged();
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

			buttonsPanelLeft.add(categories);
		}

		JPanel topPanel = new JPanel(new BorderLayout());
		topPanel.add(buttonsPanelLeft, BorderLayout.WEST);
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

	private JTable getCategoryTable(Category category) {
		// if (categoryTables.containsKey(category))
		// 	return categoryTables.get(category);

		SortableTableModel tableModel = (SortableTableModel)category.getTableModel();
		// if (tableModel == null)
		// 	tableModel = new DefaultCategoryTableModel(category);

		JTable categoryTable = new SimpleTable(manager, tableModel);
		categoryTables.put(category, categoryTable);
		return categoryTable;
	}
}
