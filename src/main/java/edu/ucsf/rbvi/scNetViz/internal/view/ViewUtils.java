package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Map;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.table.TableModel;

import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.tasks.FileCategoryTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteLeidenTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteLouvainTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteTSNETask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteUMAPTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.RemoteGraphTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.tSNETask;
import edu.ucsf.rbvi.scNetViz.internal.utils.CyPlotUtils;

public class ViewUtils {
	static int FONT_SIZE = 10;

	public static JButton addButton(ActionListener listener, String label, String command) {
		JButton button = new JButton(label);
		button.setFont(new Font("SansSerif", Font.PLAIN, FONT_SIZE));
		button.setAlignmentX(Component.CENTER_ALIGNMENT);
		button.setActionCommand(command);
		button.addActionListener(listener);
		return button;
	}

	public static JRadioButton addRadioButton(ActionListener listener, ButtonGroup group, 
	                                          String label, String command, boolean set) {
		JRadioButton button = new JRadioButton(label);
		button.setFont(new Font("SansSerif", Font.PLAIN, FONT_SIZE));
		button.setActionCommand(command);
		button.addActionListener(listener);
		button.setSelected(set);
		group.add(button);
		return button;
	}

	public static JTextField addLabeledField(JPanel container, String label, String defaultValue) {
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));
		// panel.add(Box.createHorizontalGlue());
		panel.add(Box.createRigidArea(new Dimension(15,0)));

		JLabel jlabel = new JLabel(label);
		jlabel.setFont(new Font("SansSerif", Font.BOLD, FONT_SIZE));
		jlabel.setMaximumSize(new Dimension(80,35));
		panel.add(jlabel);

		JTextField field = new JTextField(defaultValue);
		field.setFont(new Font("SansSerif", Font.PLAIN, FONT_SIZE));
		field.setMinimumSize(new Dimension(50,25));
		field.setPreferredSize(new Dimension(50,25));
		field.setMaximumSize(new Dimension(50,35));
		panel.add(field);
		panel.add(Box.createRigidArea(new Dimension(15,0)));
		container.add(panel);
		return field;
	}

	public static JPanel addJCheckBox(ActionListener listener, String label, 
	                                  String command, boolean selected) {
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));
		// panel.add(Box.createHorizontalGlue());
		panel.add(Box.createRigidArea(new Dimension(15,0)));

		JCheckBox checkBox = new JCheckBox(label, selected);
		checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
		checkBox.setFont(new Font("SansSerif", Font.BOLD, FONT_SIZE));
		checkBox.setActionCommand(command);
		checkBox.addActionListener(listener);
		panel.add(checkBox);
		panel.add(Box.createRigidArea(new Dimension(15,0)));
		return panel;
	}

	public static void showPlot(final ScNVManager manager,
                              final Experiment experiment, 
                              final Category category,
                              final int catRow,
                              final int geneRow) {
    String accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
    String title = experiment.getPlotType()+" Plot of "+accession;
		if (geneRow >= 0) {
		  TableModel tableModel = experiment.getTableModel();
			title = accession+" Expression for "+tableModel.getValueAt(geneRow, 0);
		} else if (category != null && catRow >= 0) {
			title = experiment.getPlotType()+" Plot for "+accession+" Category "+category.toString();
			List<String> rowLabels = category.getMatrix().getRowLabels();
			title += " ("+rowLabels.get(catRow)+")";
    }
		ViewUtils.showtSNE(manager, experiment, category, catRow, geneRow, title);
	}


	public static void showtSNE(final ScNVManager manager, final Experiment exp, 
	                            final Category category, final int catRow, 
	                            final int geneRow, final String title) {
		double[][] tSNEresults = exp.getTSNE();
		String type = exp.getPlotType();

		if (tSNEresults == null) {
			return;
		}

		String accession = (String)exp.getMetadata().get(Metadata.ACCESSION);
		String ttl = title;
		if (ttl == null)
			ttl = type+" Plot for "+accession;

		if (category == null) {
			// See if a gene is selected and provide a color trace if it is
			String names = "{\"trace\": "+CyPlotUtils.listToJSON(exp.getMatrix().getColLabels())+"}";
			String xValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 0)+"}";
			String yValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 1)+"}";
			String zValues = null;
			if (geneRow >= 0) {
				zValues = "{\"trace\": "+
								CyPlotUtils.valuesToJSON((DoubleMatrix)exp.getMatrix(), geneRow, true)+"}";

			}

			CyPlotUtils.createScatterPlot(manager, names, xValues, yValues, zValues, 
			                              ttl, type+" 1", type+" 2", accession);
		} else {
			String names;
			String xValues;
			String yValues;
			if (category != null && catRow >= 0) {
				Map<Object, List<Integer>> catMap = category.getCatMap(catRow);
				// Reformat the catmap so we have reasonable labels
				Map<Object, List<Integer>> newMap = new LinkedHashMap<>();

				List sortedKeys = new ArrayList(catMap.keySet());

				// Avoid cast errors
				if (sortedKeys.contains("unused"))
					sortedKeys.remove("unused");

				/*
				for (Object key: catMap.keySet()) {
					System.out.println("Key "+key.toString()+" is type "+key.getClass().getName());
				}
				*/
				Collections.sort(sortedKeys);

				for (Object key: sortedKeys) {
					newMap.put(category.mkLabel(key), catMap.get(key));
				}

				// Sort the keys now.  We need to do this now because clusters are integers and this
				// will automatically do a numeric sort.  If we waited until after we made them labels, it
				// would wind up as an alphabetical sort
				// Collections.sort(sortedKeys);

				// Now make the keys labels
				// List<String> sortedLabels = new ArrayList<String>(sortedKeys.size());
				// for (Object key: sortedKeys) {
				// 	sortedLabels.add(category.mkLabel(key));
				// }
				names = CyPlotUtils.listToMap(newMap, exp.getMatrix().getColLabels());
				xValues = CyPlotUtils.coordsToMap(newMap, tSNEresults, 0);
				yValues = CyPlotUtils.coordsToMap(newMap, tSNEresults, 1);
			} else {
				names = "{\"trace\": "+CyPlotUtils.listToJSON(exp.getMatrix().getColLabels())+"}";
				xValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 0)+"}";
				yValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 1)+"}";
			}

			CyPlotUtils.createScatterPlot(manager, names, xValues, yValues, null, 
			                              title, type+" 1", type+" 2", accession);
		}
	}

	public static JComboBox<String> createPlotMenu(final ScNVManager manager, 
	                                               final Experiment experiment, 
	                                               final TaskObserver observer) {
		Map<String, Task> map = new LinkedHashMap<>();
		String accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
		map.put("t-SNE (local)", new tSNETask(experiment));
		map.put("UMAP", new RemoteUMAPTask(manager, accession));
		map.put("Graph layout", new RemoteGraphTask(manager, accession));
		map.put("t-SNE (on server)", new RemoteTSNETask(manager, accession));
		return new PullDownMenu(manager, "New Cell Plot", map, observer);
	}

	public static JComboBox<String> createCategoryMenu(final ScNVManager manager, 
	                                                   final Experiment experiment) {
		String accession = experiment.getMetadata().get(Metadata.ACCESSION).toString();
		Map<String, Task> map = new LinkedHashMap<>();
		map.put("Import from file...",
		        new FileCategoryTask(manager, (FileSource)manager.getSource("file"), experiment));
		map.put("Louvain clustering...",
		        new RemoteLouvainTask(manager, accession));
		map.put("Leiden clustering...",
		        new RemoteLeidenTask(manager, accession));
		return new PullDownMenu(manager, "Add Category", map, null);
	}

}
