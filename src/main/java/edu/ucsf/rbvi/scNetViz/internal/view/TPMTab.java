package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.ScrollPaneConstants;
import javax.swing.table.TableModel;

import java.util.HashMap;
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
import edu.ucsf.rbvi.scNetViz.internal.tasks.ExportCSVTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.tSNETask;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;

public class TPMTab extends JPanel implements TaskObserver {
	final ScNVManager manager;
	final Experiment experiment;
	final TPMTab thisComponent;
	final ExperimentFrame expFrame;
	JTable experimentTable;

	public TPMTab(final ScNVManager manager, final Experiment experiment, final ExperimentFrame expFrame) {
		this.manager = manager;
		this.experiment = experiment;

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
			int index = rowLabels.indexOf(gene);
			experimentTable.getSelectionModel().addSelectionInterval(index, index);
		}
	}

	@Override
	public void allFinished(FinishStatus status) {
	}

	@Override
	public void taskFinished(ObservableTask obsTask) {
		if (obsTask instanceof FileCategoryTask) {
			String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
			expFrame.addCategoriesContent(accession+": Categories Tab", new CategoriesTab(manager, experiment, expFrame));
		} else if (obsTask instanceof tSNETask) {
			double[][] tSNEResults = ((tSNETask)obsTask).getResults();
			experiment.setTSNE(tSNEResults);
			int geneRow = experimentTable.getSelectedRow();
			String title = null;
			if (geneRow >= 0) {
				String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
				geneRow = experimentTable.convertRowIndexToModel(geneRow);
				title = accession+" Gene "+experimentTable.getModel().getValueAt(geneRow, 0);
			}
			ViewUtils.showtSNE(manager, experiment, null, -1, geneRow, title);
		}
	}
	
	private void init() {

		// JLabel experimentLabel = new ExperimentLabel(experiment);

		JPanel buttonsPanelRight = new JPanel();
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

		{
			JButton tsne = new JButton("View tSNE");
			tsne.setFont(new Font("SansSerif", Font.PLAIN, 10));
      tsne.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					String title = null;
					int geneRow = experimentTable.getSelectedRow();
					if (geneRow >= 0) {
						String accession = (String)experiment.getMetadata().get(Metadata.ACCESSION);
						geneRow = experimentTable.convertRowIndexToModel(geneRow);
						title = accession+" Expression for "+experimentTable.getModel().getValueAt(geneRow, 0);
					}
					ViewUtils.showtSNE(manager, experiment, null, -1, geneRow, title);
				}
			});
			buttonsPanelRight.add(tsne);
		}

		{
			JButton tsne = new JButton("(Re)calculate tSNE");
			tsne.setFont(new Font("SansSerif", Font.PLAIN, 10));
      tsne.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					Task tSNETask = new tSNETask((DoubleMatrix)experiment.getMatrix());
					manager.executeTasks(new TaskIterator(tSNETask), thisComponent);
				}
			});
			buttonsPanelRight.add(tsne);
		}
		
		{
			JButton importCategory = new JButton("Import Category");
			importCategory.setFont(new Font("SansSerif", Font.PLAIN, 10));
      importCategory.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					// We need to use the File importer for this
					TaskFactory importCategory = 
									new FileCategoryTaskFactory(manager, (FileSource)manager.getSource("file"), experiment);
					manager.executeTasks(importCategory, thisComponent);
				}
			});
			buttonsPanelRight.add(importCategory);
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

		JScrollPane scrollPane = new JScrollPane(experimentTable);
		scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		this.add(scrollPane, BorderLayout.CENTER);
		this.revalidate();
		this.repaint();
	}

}
