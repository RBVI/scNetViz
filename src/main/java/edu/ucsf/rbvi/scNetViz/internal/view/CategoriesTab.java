package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ScrollPaneConstants;
import javax.swing.table.TableModel;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class CategoriesTab extends JPanel {
	final ScNVManager manager;
	final Experiment experiment;
	final CategoriesTab thisComponent;
	final List<Category> categoriesList;
	final List<String> categoriesNameList;
	final Map<Category, JTable> categoryTables;
	JScrollPane categoryPane;

	public CategoriesTab(final ScNVManager manager, final Experiment experiment) {
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
		init();
	}

	private void init() {

		JPanel buttonsPanelRight = new JPanel(new GridLayout(3, 2));
		
		{
			JLabel lbl = new JLabel("");
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
				}
			});

			buttonsPanelRight.add(lbl);
			buttonsPanelRight.add(export);
			buttonsPanelRight.add(diffExp);
		}
		
		JPanel buttonsPanelLeft = new JPanel(new GridLayout(3,1));
		{
			JLabel experimentLabel = new ExperimentLabel(experiment);
			JLabel lbl = new JLabel("Available categories");
			lbl.setFont(new Font("SansSerif", Font.PLAIN, 10));
			JComboBox categories = new JComboBox(categoriesNameList.toArray());
      categories.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					System.out.println("Switch table to: "+categoriesList.get(categories.getSelectedIndex()));
					JTable table = getCategoryTable(categoriesList.get(categories.getSelectedIndex()));
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
			buttonsPanelLeft.add(experimentLabel);
			buttonsPanelLeft.add(lbl);
			buttonsPanelLeft.add(categories);
		}

		JPanel topPanel = new JPanel(new BorderLayout());
		topPanel.add(buttonsPanelLeft, BorderLayout.WEST);
		topPanel.add(buttonsPanelRight, BorderLayout.EAST);
		this.add(topPanel, BorderLayout.NORTH);

		JTable categoryTable = getCategoryTable(categoriesList.get(0));

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
