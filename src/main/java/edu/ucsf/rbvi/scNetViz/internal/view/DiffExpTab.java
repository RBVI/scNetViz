package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.table.TableModel;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class DiffExpTab extends JPanel {
	final ScNVManager manager;
	final Experiment experiment;
	final ExperimentFrame expFrame;
	final List<Category> categories;
	final Category currentCategory;
	final DifferentialExpression diffExp;
	final DiffExpTab thisComponent;
	final Map<Category, List<String>> categoryLabelMap;

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

	private void init() {

		JPanel buttonsPanelRight = new JPanel(new GridLayout(4, 1));
		
		{
			JButton viewHeatMap = new JButton("View Heatmap");
			viewHeatMap.setFont(new Font("SansSerif", Font.PLAIN, 10));
      viewHeatMap.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
				}
			});

			JButton export = new JButton("Export CSV");
			export.setFont(new Font("SansSerif", Font.PLAIN, 10));
      export.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
				}
			});

			JButton viewViolin = new JButton("View Violin Plot");
			viewViolin.setFont(new Font("SansSerif", Font.PLAIN, 10));
      viewViolin.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
				}
			});

			buttonsPanelRight.add(new JLabel(""));
			buttonsPanelRight.add(viewHeatMap);
			buttonsPanelRight.add(viewViolin);
			buttonsPanelRight.add(export);
		}
		
		JPanel buttonsPanelLeft = new JPanel();
		buttonsPanelLeft.setLayout(new BoxLayout(buttonsPanelLeft, BoxLayout.PAGE_AXIS));
		{
			{
				JLabel experimentLabel = new ExperimentLabel(experiment);
				experimentLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
				buttonsPanelLeft.add(experimentLabel);
				buttonsPanelLeft.add(Box.createRigidArea(new Dimension(0, 10)));
			}

			{
				JLabel lbl = new JLabel("Comparison:");
				lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				lbl.setAlignmentX(Component.LEFT_ALIGNMENT);
				buttonsPanelLeft.add(lbl);
			}

			{
				List<String> labels = new ArrayList<>();
				for (List<String> lbl: categoryLabelMap.values())
					labels.addAll(lbl);

				int selectedRow = currentCategory.getSelectedRow();
				String selectedLabel = categoryLabelMap.get(currentCategory).get(selectedRow);

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
						// Reset the selected row
						// Recalculate
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
			{
				JLabel lbl = new JLabel("Network analysis:");
				lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				lbl.setAlignmentX(Component.LEFT_ALIGNMENT);
				centerPanel.add(lbl);
			}

			JPanel settingsPanel = new JPanel();
			settingsPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
			settingsPanel.setLayout(new BoxLayout(settingsPanel, BoxLayout.LINE_AXIS));
			{
				settingsPanel.add(Box.createRigidArea(new Dimension(5,0)));
				JLabel pValueLbl = new JLabel("P-value");
				pValueLbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				pValueLbl.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(pValueLbl);
			}

			{
				JTextField pValue = new JTextField("0.05");
				pValue.setFont(new Font("SansSerif", Font.PLAIN, 10));
				pValue.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(pValue);
				settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));
			}

			{
				JLabel log2FCLbl = new JLabel("Log2FC");
				log2FCLbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				log2FCLbl.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(log2FCLbl);
			}

			{
				JTextField log2FC = new JTextField("1.0");
				log2FC.setFont(new Font("SansSerif", Font.PLAIN, 10));
				log2FC.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(log2FC);
				settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));
			}

			{
				JLabel topGenesLbl = new JLabel("Top n genes");
				topGenesLbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				topGenesLbl.setMaximumSize(new Dimension(70,35));
				settingsPanel.add(topGenesLbl);
			}

			{
				JTextField topGenes = new JTextField("25");
				topGenes.setFont(new Font("SansSerif", Font.PLAIN, 10));
				topGenes.setMaximumSize(new Dimension(50,35));
				settingsPanel.add(topGenes);
				settingsPanel.add(Box.createRigidArea(new Dimension(15,0)));
			}

			{
				JButton createNetwork = new JButton("Create Network");
				createNetwork.setFont(new Font("SansSerif", Font.BOLD, 10));
				createNetwork.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {

					}
				});
				settingsPanel.add(createNetwork);
			}

			centerPanel.add(settingsPanel);

		}


		JPanel topPanel = new JPanel(new BorderLayout());
		topPanel.add(buttonsPanelLeft, BorderLayout.WEST);
		topPanel.add(centerPanel, BorderLayout.CENTER);
		topPanel.add(buttonsPanelRight, BorderLayout.EAST);
		this.add(topPanel, BorderLayout.NORTH);

		SortableTableModel tableModel = diffExp.getTableModel();

		JTable diffExpTable = new SimpleTable(manager, tableModel);

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
}
