package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableModel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAExperiment;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTask;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.tasks.AbstractEmbeddingTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ExportCSVTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteTSNETask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteUMAPTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteGraphTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.tSNETask;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;
import edu.ucsf.rbvi.scNetViz.internal.utils.ModelUtils;

public class TPMTab extends JPanel implements TaskObserver {
	final ScNVManager manager;
	final Experiment experiment;
	final TPMTab thisComponent;
	final ExperimentFrame expFrame;
	final String accession;
	JTable experimentTable;
	public JButton cellPlotButton;

	public TPMTab(final ScNVManager manager, final Experiment experiment, final ExperimentFrame expFrame) {
		this.manager = manager;
		this.experiment = experiment;
		this.accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();

		this.setLayout(new BorderLayout());
		thisComponent = this;	// Access to inner classes
		this.expFrame = expFrame;
		init();
	}

	public void selectGenes(List<String> geneList) {
		// Clear the selection list
		experimentTable.getSelectionModel().clearSelection();
		// Get the unsorted row labels
		List<String> rowLabels = experiment.getMatrix().getRowLabels();
		for (String gene: geneList) {
			System.out.println("Selecting gene: "+gene);
			int index = rowLabels.indexOf(gene);
			index = experimentTable.convertRowIndexToView(index);
			experimentTable.getSelectionModel().addSelectionInterval(index, index);
			experimentTable.scrollRectToVisible(new Rectangle(experimentTable.getCellRect(index, 0, true)));
		}
		ModelUtils.selectNodes(manager, accession, geneList);
	}

	@Override
	public void allFinished(FinishStatus status) {
	}

	@Override
	public void taskFinished(ObservableTask obsTask) {
		if (obsTask instanceof FileCategoryTask) {
			expFrame.addCategoriesContent(accession+": Categories Tab", new CategoriesTab(manager, experiment, expFrame));
		} else if (obsTask instanceof AbstractEmbeddingTask) {
			double[][] embedding = ((AbstractEmbeddingTask)obsTask).getResults();
			if (embedding == null)
				return;
			experiment.setTSNE(embedding);
			showPlot();
			cellPlotButton.setEnabled(true);
			cellPlotButton.setText("View "+experiment.getPlotType());
			CategoriesTab cTab = expFrame.getCategoriesTab();
			if (cTab != null) {
				cTab.cellPlotButton.setEnabled(true);
				cTab.cellPlotButton.setText("View "+experiment.getPlotType());
			}
			if (manager.getCytoPanel() != null) {
				manager.getCytoPanel().updatePlotMenu();
			}
		}
		expFrame.toFront();
	}
	
	private void init() {

		// JLabel experimentLabel = new ExperimentLabel(experiment);

		JPanel buttonsPanelRight = new JPanel();

		{
			buttonsPanelRight.add(ViewUtils.createPlotMenu(manager, experiment, thisComponent));
		}

		{
			// TODO: convert to pull down with "Plots", "Cell Plot", "Gene Plot"
			cellPlotButton = new JButton("View Cell Plot");
			cellPlotButton.setEnabled(false);
			cellPlotButton.setFont(new Font("SansSerif", Font.PLAIN, 10));
      cellPlotButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					showPlot();
				}
			});
			buttonsPanelRight.add(cellPlotButton);
		}
		
		{
			buttonsPanelRight.add(ViewUtils.createCategoryMenu(manager, experiment));
		}

		{
			JButton export = new JButton("Export CSV");
			export.setFont(new Font("SansSerif", Font.PLAIN, 10));
      export.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					ExportCSVTask task = new ExportCSVTask(manager, experiment.getMatrix());
					manager.executeTasks(new TaskIterator(task));
				}
			});
			buttonsPanelRight.add(export);
		}
		
		JPanel topPanel = new JPanel(new BorderLayout());
		//topPanel.add(experimentLabel, BorderLayout.WEST);
		topPanel.add(buttonsPanelRight, BorderLayout.EAST);
		this.add(topPanel, BorderLayout.NORTH);

		TableModel tableModel = experiment.getTableModel();
		if (tableModel == null)
			tableModel = new DefaultExperimentTableModel(experiment);

		experimentTable = new ExperimentTable(manager, tableModel);
		experimentTable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		experimentTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent event) {
				int[] rows = experimentTable.getSelectedRows();
				if (rows.length == 0) return;
				List<String> geneList = new ArrayList<>();
				for (int row: rows) {
					geneList.add(experimentTable.getValueAt(row, 0).toString());
				}
				ModelUtils.selectNodes(manager, accession, geneList);
			}
		});

		JScrollPane scrollPane = new JScrollPane(experimentTable);
		scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		this.add(scrollPane, BorderLayout.CENTER);
		this.revalidate();
		this.repaint();
	}

	private void showPlot() {
		String title = null;
		int geneRow = experimentTable.getSelectedRow();
		if (geneRow >= 0) {
			geneRow = experimentTable.convertRowIndexToModel(geneRow);
			title = accession+" Expression for "+experimentTable.getModel().getValueAt(geneRow, 0);
		}
		ViewUtils.showtSNE(manager, experiment, null, -1, geneRow, title);
	}

}
