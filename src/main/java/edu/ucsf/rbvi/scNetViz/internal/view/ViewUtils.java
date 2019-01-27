package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.HashMap;
import java.util.List;
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
		panel.add(Box.createHorizontalGlue());

		JLabel jlabel = new JLabel(label);
		jlabel.setFont(new Font("SansSerif", Font.BOLD, FONT_SIZE));
		jlabel.setMaximumSize(new Dimension(100,35));
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
		panel.add(Box.createHorizontalGlue());

		JCheckBox checkBox = new JCheckBox(label, selected);
		checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
		checkBox.setFont(new Font("SansSerif", Font.BOLD, FONT_SIZE));
		checkBox.setActionCommand(command);
		checkBox.addActionListener(listener);
		panel.add(checkBox);
		panel.add(Box.createRigidArea(new Dimension(15,0)));
		return panel;
	}

	public static void showtSNE(final ScNVManager manager, final Experiment exp, 
	                            final Category category, final int catRow, 
	                            final int geneRow, final String title) {
		double[][] tSNEresults = exp.getTSNE();

		if (tSNEresults == null) {
			Task tSNETask = new tSNETask((DoubleMatrix)exp.getMatrix());
			manager.executeTasks(new TaskIterator(tSNETask), new TaskObserver() {
				@Override
				public void allFinished(FinishStatus status) {
				}

				@Override
				public void taskFinished(ObservableTask obsTask) {
					if (obsTask instanceof tSNETask) {
						double[][] tSNEResults = ((tSNETask)obsTask).getResults();
						exp.setTSNE(tSNEResults);
						showtSNE(manager, exp, category, catRow, geneRow, title);
					}
				}
			});
			return;
		}

		String accession = (String)exp.getMetadata().get(Metadata.ACCESSION);
		String ttl = title;
		if (ttl == null)
			ttl = "tSNE Plot for "+accession;

		if (category == null) {
			// See if a gene is selected and provide a color trace if it is
			String names = "{\"trace\": "+CyPlotUtils.listToJSON(exp.getMatrix().getColLabels())+"}";
			String xValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 0)+"}";
			String yValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 1)+"}";
			String zValues = null;
			if (geneRow >= 0)
				zValues = "{\"trace\": "+
								CyPlotUtils.valuesToJSON((DoubleMatrix)exp.getMatrix(), geneRow)+"}";

			CyPlotUtils.createScatterPlot(manager, names, xValues, yValues, zValues, 
			                              ttl, "t-SNE 1", "t-SNE 2", accession);
		} else {
			String names;
			String xValues;
			String yValues;
			if (category != null && catRow >= 0) {
				Map<Object, List<Integer>> catMap = category.getCatMap(catRow);
				// Reformat the catmap so we have reasonable labels
				Map<Object, List<Integer>> newMap = new HashMap<>();
				for (Object key: catMap.keySet()) {
					if (key.toString().equals("unused"))
						continue;
					newMap.put(category.mkLabel(key), catMap.get(key));
				}
				names = CyPlotUtils.listToMap(newMap, exp.getMatrix().getColLabels());
				xValues = CyPlotUtils.coordsToMap(newMap, tSNEresults, 0);
				yValues = CyPlotUtils.coordsToMap(newMap, tSNEresults, 1);
			} else {
				names = "{\"trace\": "+CyPlotUtils.listToJSON(exp.getMatrix().getColLabels())+"}";
				xValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 0)+"}";
				yValues = "{\"trace\": "+CyPlotUtils.coordinatesToJSON(tSNEresults, 1)+"}";
			}

			CyPlotUtils.createScatterPlot(manager, names, xValues, yValues, null, 
			                              title, "t-SNE 1", "t-SNE 2", accession);
		}
	}

}
