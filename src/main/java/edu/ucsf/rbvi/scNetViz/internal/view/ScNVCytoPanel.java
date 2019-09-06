package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.events.SetCurrentNetworkListener;
import org.cytoscape.application.events.SetCurrentNetworkEvent;
import org.cytoscape.application.swing.CytoPanelComponent2;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.CyIdentifiable;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.CyTableUtil;
import org.cytoscape.model.events.RowSetRecord;
import org.cytoscape.model.events.RowsSetEvent;
import org.cytoscape.model.events.RowsSetListener;
import org.cytoscape.util.swing.IconManager;
import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVSettings.SETTING;
import edu.ucsf.rbvi.scNetViz.internal.tasks.CreateNetworkTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.HeatMapTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ViolinDiffExpTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ShowExperimentTableTask;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.LogUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

/**
 * Displays controls for scNetViz
 * @author Scooter Morris
 *
 */
public class ScNVCytoPanel extends JPanel 
                           implements CytoPanelComponent2, RowsSetListener, ActionListener,
                                      SetCurrentNetworkListener {

	final ScNVManager manager;
	final CyApplicationManager appManager;
	final Font iconFont;
	private List<Category> categories;
	private Map<Category, List<String>> categoryLabelMap;
	private String accession;
	private String currentCategoryColumn;

	private CyNetwork network;
	private Experiment experiment;

	private JTextField log2FC;
	private JTextField pValue;
	private JTextField topNGenes;
	private JTextField maxGenes;
	private JTextField FDRCutoff;
	private ExperimentLabel experimentLabel;

	private String enrichmentType = "entireNetwork";

	private boolean positiveOnly = false;

	private boolean tSNEShown = false;

	private PullDownMenu plotMenu;

	public ScNVCytoPanel(final ScNVManager manager, final Experiment experiment) {
		this.manager = manager;
		this.appManager = manager.getService(CyApplicationManager.class);
		this.network = appManager.getCurrentNetwork();
		this.experiment = experiment;
		if (this.experiment == null) {
			this.experiment = ModelUtils.getExperimentFromNetwork(manager, network);
		}

		if (this.experiment == null)
			throw new RuntimeException("Experiment is null!");

		accession = this.experiment.getMetadata().get(Metadata.ACCESSION).toString();
		categories = this.experiment.getCategories();

		manager.setCytoPanel(this);

		categoryLabelMap = new HashMap<>();
		for (Category cat: categories) {
			categoryLabelMap.put(cat, cat.getMatrix().getRowLabels());
		}
		currentCategoryColumn = ModelUtils.getName(network, network);

		IconManager iconManager = manager.getService(IconManager.class);
		iconFont = iconManager.getIconFont(17.0f);
		setLayout(new BorderLayout());
		add(createLabelPanel(), BorderLayout.NORTH);
		JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.PAGE_AXIS));
		mainPanel.add(createViewPanel());
		mainPanel.add(createComparisonPanel());
		mainPanel.add(createEnrichmentPanel());
		// TODO: add a heat map for currently selected nodes
		// mainPanel.add(createHeatMapPanel());
		JScrollPane scrollPane = new JScrollPane(mainPanel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
		                                         JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		add(scrollPane, BorderLayout.CENTER);
		revalidate();
		repaint();
	}

	public String getIdentifier() {
		return "edu.ucsf.rbvi.scNetViz.ResultsPanel";
	}

	@Override
	public void handleEvent(RowsSetEvent rse) {
		// TODO: implement when we implement ResultsPanel heatmap
	}

	@Override
	public void handleEvent(SetCurrentNetworkEvent scne) {
		network = scne.getNetwork();
		experiment = ModelUtils.getExperimentFromNetwork(manager, network);
		if (experiment == null) return;
		setExperiment(experiment);
	}

	public Experiment getExperiment() {
		return experiment;
	}

	public void setExperiment(Experiment experiment) {
		this.experiment = experiment;
		accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
		categories = experiment.getCategories();

		categoryLabelMap = new HashMap<>();
		for (Category cat: categories) {
			categoryLabelMap.put(cat, cat.getMatrix().getRowLabels());
		}
		experimentLabel.updateText(experiment);
		currentCategoryColumn = ModelUtils.getName(network, network);
	}

	@Override
	public Component getComponent() {
		// TODO Auto-generated method stub
		return this;
	}

	@Override
	public CytoPanelName getCytoPanelName() {
		// TODO Auto-generated method stub
		return CytoPanelName.EAST;
	}

	@Override
	public Icon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getTitle() {
		return "ScNetViz";
	}

	public void updatePlotMenu() {
		plotMenu.updateMenu(getPlotMap());
	}

	private Map<String, Task> getPlotMap() {
		Map<String, Task> menu = new LinkedHashMap<>();
		String plotType = experiment.getPlotType();
		if (plotType != null) {
			menu.put(plotType, new CellPlot());
		}
		menu.put("Heatmap", new HeatMapWrapperTask());
		menu.put("Violin", new ViolinWrapperTask());
		return menu;
	}


	private JPanel createLabelPanel() {
		JPanel labelPanel = new JPanel();
		labelPanel.setLayout(new BoxLayout(labelPanel, BoxLayout.PAGE_AXIS));
		// labelPanel.add(new JSeparator(SwingConstants.HORIZONTAL));
		{
			JPanel pnl = new JPanel();
			pnl.setLayout(new BoxLayout(pnl, BoxLayout.LINE_AXIS));
			experimentLabel = new ExperimentLabel(experiment, pnl.getBackground(), 150);
			pnl.add(Box.createRigidArea(new Dimension(2,0)));
			pnl.add(experimentLabel);
			labelPanel.add(pnl);
		}
		{
			JPanel pnl = new JPanel();
			pnl.setLayout(new BoxLayout(pnl, BoxLayout.LINE_AXIS));
			pnl.add(Box.createHorizontalGlue());
			pnl.add(new HelpButton(manager, null));
			labelPanel.add(pnl);
		}
		return labelPanel;
	}

	private JPanel createViewPanel() {
		JPanel viewPanel = new JPanel();
		Dimension size = new Dimension(200,25);
		viewPanel.setLayout(new BoxLayout(viewPanel, BoxLayout.PAGE_AXIS));
		viewPanel.add(Box.createHorizontalGlue());
		Map<String, Task> menuMap = new LinkedHashMap<>();
		menuMap.put("TPM Table", new ShowExperimentTableTask(manager, experiment, "tpmTable"));
		menuMap.put("Category Table", new ShowExperimentTableTask(manager, experiment, "catTable"));
		menuMap.put("DE Table", new ShowExperimentTableTask(manager, experiment, "DETable"));
		PullDownMenu tableMenu = new PullDownMenu(manager, "Tables", menuMap, null);
		tableMenu.setPreferredSize(size);
		tableMenu.setMaximumSize(size);
		tableMenu.setMinimumSize(size);
		viewPanel.add(tableMenu);

		plotMenu = new PullDownMenu(manager, "Plots", getPlotMap(), null);
		plotMenu.setPreferredSize(size);
		plotMenu.setMaximumSize(size);
		plotMenu.setMinimumSize(size);
		viewPanel.add(plotMenu);
		// viewPanel.add(Box.createHorizontalGlue());
		viewPanel.add(Box.createVerticalGlue());
		return new CollapsablePanel(iconFont, "View", viewPanel, false);
	}

	private JPanel createComparisonPanel() {
		JPanel comparePanel = new JPanel();
		comparePanel.setLayout(new BoxLayout(comparePanel, BoxLayout.PAGE_AXIS));
		{
			{
				JPanel panel = new JPanel();
				panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));

				JLabel lbl = new JLabel("Comparison:");
				lbl.setFont(new Font("SansSerif", Font.BOLD, 10));
				panel.add(Box.createRigidArea(new Dimension(5,5)));
				panel.add(lbl);
				panel.add(Box.createHorizontalGlue());
				comparePanel.add(panel);
			}

			{
				JPanel panel = new JPanel();
				panel.setLayout(new BoxLayout(panel, BoxLayout.PAGE_AXIS));
				// JPanel panel = new JPanel();
				// panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));
				// panel.add(Box.createRigidArea(new Dimension(5,0)));
				panel.add(Box.createHorizontalGlue());

				List<String> labels = new ArrayList<>();
				for (List<String> lbl: categoryLabelMap.values())
					labels.addAll(lbl);

				// int selectedRow = currentCategory.getSelectedRow();
				// String selectedLabel = categoryLabelMap.get(currentCategory).get(selectedRow);

				JComboBox<String> categoryBox = 
					new JComboBox<String>(labels.toArray(new String[1]));

				String categoryRow = ModelUtils.getCategoryRowFromNetwork(manager, network);
				System.out.println("CurrentCategoryColumn="+categoryRow);

				categoryBox.setSelectedItem(categoryRow);

			 	Dimension size = new Dimension(200,25);
				categoryBox.setPreferredSize(size);
				categoryBox.setMaximumSize(new Dimension(250,25));
				categoryBox.setSize(size);
				categoryBox.setFont(new Font("SansSerif", Font.PLAIN, 10));
				categoryBox.setAlignmentX(Component.LEFT_ALIGNMENT);
				panel.add(categoryBox);
				panel.add(Box.createHorizontalGlue());
				// panel.add(Box.createHorizontalGlue());
				panel.add(Box.createRigidArea(new Dimension(5,0)));

				JComboBox comparison = new JComboBox(getComparisons());
				comparison.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
					}
				});
				comparison.setPreferredSize(size);
				comparison.setMaximumSize(new Dimension(250,25));
				comparison.setSize(size);
				comparison.setFont(new Font("SansSerif", Font.PLAIN, 10));
				comparison.setAlignmentX(Component.LEFT_ALIGNMENT);
				panel.add(comparison);
				comparePanel.add(panel);
				// panel.add(Box.createHorizontalGlue());
			}

			pValue = ViewUtils.addLabeledField(comparePanel, "pValue:", 
			                                   manager.getSetting(SETTING.NET_PV_CUTOFF));
			log2FC = ViewUtils.addLabeledField(comparePanel, "log2FC:", 
			                                   manager.getSetting(SETTING.NET_FC_CUTOFF));
			topNGenes = ViewUtils.addLabeledField(comparePanel, "Top n genes:", 
			                                      manager.getSetting(SETTING.TOP_GENES));
			maxGenes = ViewUtils.addLabeledField(comparePanel, "Max genes:", 
			                                     manager.getSetting(SETTING.MAX_GENES));
			comparePanel.add(ViewUtils.addJCheckBox(this, "Positive only", "positiveOnlyComp", positiveOnly));
			comparePanel.add(ViewUtils.addButton(this, "Create Networks", "createNetworks"));
			// comparePanel.add(Box.createRigidArea(new Dimension(0, 15)));
			// comparePanel.add(ViewUtils.addButton(this, "View Heatmap", "viewHeatmap"));
			// comparePanel.add(ViewUtils.addButton(this, "View Violin Plot", "viewViolin"));
		}

		return new CollapsablePanel(iconFont, "Reanalyze", comparePanel, false);
	}

	private JPanel createEnrichmentPanel() {
		JPanel enrichmentPanel = new JPanel();
		enrichmentPanel.setLayout(new BoxLayout(enrichmentPanel, BoxLayout.PAGE_AXIS));
		{
			ButtonGroup group = new ButtonGroup();
			enrichmentPanel.add(ViewUtils.addRadioButton(this, group,
			                                             "Entire network", ModelUtils.ENTIRE_NETWORK, true));
			enrichmentPanel.add(ViewUtils.addRadioButton(this, group,
			                                             "Positive only", ModelUtils.POSITIVE_ONLY, false));
			enrichmentPanel.add(ViewUtils.addRadioButton(this, group,
			                                             "Negative only", ModelUtils.NEGATIVE_ONLY, false));
			enrichmentPanel.add(ViewUtils.addRadioButton(this, group,
			                                             "Selected only", ModelUtils.SELECTED_ONLY, false));
			FDRCutoff = ViewUtils.addLabeledField(enrichmentPanel, "FDR cutoff:", "0.05");

			enrichmentPanel.add(Box.createRigidArea(new Dimension(0, 5)));
			{
				JPanel panel = new JPanel();
				panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));
				// panel.add(Box.createRigidArea(new Dimension(5,0)));
				panel.add(Box.createHorizontalGlue());
				panel.add(ViewUtils.addButton(this, "Retrieve Table", "getEnrichment"));
				panel.add(Box.createHorizontalGlue());
				enrichmentPanel.add(panel);
			}
			enrichmentPanel.add(Box.createRigidArea(new Dimension(0, 5)));

		}

		return new CollapsablePanel(iconFont, "Get Enrichment", enrichmentPanel, false);
	}

	public void actionPerformed(ActionEvent event) {
		String command = event.getActionCommand();
		switch (command) {
			case "positiveOnlyComp":
				{
					JCheckBox po = (JCheckBox)event.getSource();
					positiveOnly = po.isSelected();
					break;
				}

			case "entireNetwork":
			case "positiveOnly":
			case "negativeOnly":
			case "selectedOnly":
				enrichmentType = command;
				break;

			case "getEnrichment":
				{
					// Perform enrichment
					double fdrCutoff = Double.parseDouble(FDRCutoff.getText());

					List<CyNode> currentSelection = null;

					String selectedOnly = "false";
					if (enrichmentType != ModelUtils.ENTIRE_NETWORK) {
						if (enrichmentType != ModelUtils.SELECTED_ONLY)
							currentSelection = ModelUtils.selectNodes(manager, network, enrichmentType);
						selectedOnly = "true";
					}
					// Call string get enrichment command
					Map<String, Object> args = new HashMap<>();

					args.clear();
					// Hide the enrichment panel (just in case)
					manager.executeCommand("string", "hide enrichment", args, null, true);

					args.clear();
					args.put("cutoff", String.valueOf(fdrCutoff));
					args.put("selectedNodesOnly", selectedOnly);
					manager.executeCommand("string", "retrieve enrichment", args, null, true);

					args.clear();
					// Show the enrichment panel
					manager.executeCommand("string", "show enrichment", args, null, true);

					args.clear();
					args.put("removeOverlapping", "true");
					// Filter the enrichment data
					manager.executeCommand("string", "filter enrichment", args, null, true);

					args.clear();
					// Show our charts
					manager.executeCommand("string", "show charts", args, null, true);

					args.clear();
					// Hide the glass effect
					manager.executeCommand("string", "show glass", args, null, true);

					args.clear();
					// Hide the glass effect
					manager.executeCommand("string", "hide glass", args, null, true);

					if (currentSelection != null) {
						ModelUtils.selectNodes(network, currentSelection);
					}
				}
				break;

			case "createNetworks":
				{
					DifferentialExpression diffExp = experiment.getDiffExp();
					if (diffExp == null) {
						// TODO: calculate differential expression if none exists
						LogUtils.error("No differential expression calculated");
					}
					double pV = Double.parseDouble(pValue.getText());
					double log2FCCutoff = Double.parseDouble(log2FC.getText());
					int nGenes = -1;
				 	if (topNGenes.getText() != null && !topNGenes.getText().equals(""))
						nGenes = Integer.parseInt(topNGenes.getText());
					int maxG = Integer.parseInt(maxGenes.getText());
					Task createNetworks = new CreateNetworkTask(manager, diffExp, pV, log2FCCutoff, 
					                                            nGenes, positiveOnly, maxG); 
					manager.executeTasks(new TaskIterator(createNetworks));
				}
				break;
			default:
		}
	}

	private JPanel createHeatMapPanel() {
		JPanel heatMapPanel = new JPanel();
		return new CollapsablePanel(iconFont, "HeatMap", heatMapPanel, false);
	}

	String[] getComparisons() {
		String[] comparisonTypes = {"Each vs. Others"};
		return comparisonTypes;
	}

	class CellPlot extends AbstractTask {
		public void run(TaskMonitor monitor) {
			List<CyNode> selectedNodes = ModelUtils.getSelectedNodes(network);
			int geneRow = ModelUtils.getRowFromNode(experiment, network, selectedNodes);
			String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
			String plotType = experiment.getPlotType();
			if (plotType == null)
				return; // Shouldn't happen!
			String title = plotType+" Plot for "+accession;
			if (selectedNodes != null && selectedNodes.size() > 0)
				title = accession+" Expression for "+ModelUtils.getGeneNameFromNode(network, selectedNodes.get(0));
			ViewUtils.showtSNE(manager, experiment, null, -1, geneRow, title);
		}
	}

	class HeatMapWrapperTask extends AbstractTask {
		public void run(TaskMonitor monitor) {
			// Get the category
			Category currentCategory = ModelUtils.getCategoryFromNetwork(manager, network);
			List<String> categoryNames = ModelUtils.getCategoryNamesFromNetwork(network);

			// Get the selected nodes (if any)
			List<CyNode> selectedNodes = ModelUtils.getSelectedNodes(network);
			if (selectedNodes == null || selectedNodes.size() == 0)
				selectedNodes = network.getNodeList();

			List<String> geneNames = ModelUtils.getGeneNamesFromNetwork(network, selectedNodes);
			Map<String, double[]> dataMap = new HashMap<>();
			for (String cat: categoryNames) {
				dataMap.put(cat, ModelUtils.getDataFromNetwork(network, cat+" log2FC", selectedNodes));
			}

			// Use a separate task for this since we've got some options...
			HeatMapTask task = new HeatMapTask(manager, currentCategory, geneNames, dataMap, categoryNames, 
			                                   positiveOnly, selectedNodes.size(), currentCategoryColumn);
			manager.executeTasks(new TaskIterator(task));
		}
	}

	class ViolinWrapperTask extends AbstractTask {
		public void run(TaskMonitor monitor) {
			// Get the category
			Category currentCategory = ModelUtils.getCategoryFromNetwork(manager, network);
			List<String> categoryNames = ModelUtils.getCategoryNamesFromNetwork(network);

			// Get the selected nodes (if any)
			List<CyNode> selectedNodes = ModelUtils.getSelectedNodes(network);
			if (selectedNodes == null || selectedNodes.size() == 0)
				selectedNodes = network.getNodeList();

			List<String> geneNames = ModelUtils.getGeneNamesFromNetwork(network, selectedNodes);
			Map<String, double[]> dataMap = new HashMap<>();
			for (String cat: categoryNames) {
				dataMap.put(cat, ModelUtils.getDataFromNetwork(network, cat+" log2FC", selectedNodes));
			}

			String[] dataAndNames = CyPlotUtils.mapToDataAndNames(dataMap, geneNames, categoryNames);
			String accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
			String title = experiment.getSource().toString()+" "+ accession+ " Differential Expression";
			CyPlotUtils.createViolinPlot(manager, dataAndNames[0], dataAndNames[1], 
			                             CyPlotUtils.listToCSV(categoryNames), title, "", "Log(FC)", accession, false);
		}
	}
}
