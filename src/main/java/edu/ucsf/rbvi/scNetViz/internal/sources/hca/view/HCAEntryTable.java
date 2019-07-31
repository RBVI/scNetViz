package edu.ucsf.rbvi.scNetViz.internal.sources.hca.view;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import java.util.List;
import java.util.Properties;

import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;

import org.cytoscape.application.swing.CySwingApplication;
import org.cytoscape.application.swing.CytoPanel;
import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.CyTable;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCAExperiment;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCASource;
import edu.ucsf.rbvi.scNetViz.internal.view.HeaderRenderer;

public class HCAEntryTable extends JTable {
	final ScNVManager scNVManager;
	final HCASource hcaSource;
	List<Metadata> entries;

	static String[] columnNames = {"Loaded", "Experiment", "Cells", "Comparisons", "Organisms", "Experimental Variables"};

	static Color alternateColor = new Color(234,255,234);

	public HCAEntryTable (final ScNVManager scNVManager, final HCASource hcaSource, final HCAEntryTableModel tableModel) {
		super(tableModel);
		this.scNVManager = scNVManager;
		this.hcaSource = hcaSource;

		entries = hcaSource.getMetadata();

		this.setAutoCreateRowSorter(true);
		this.setAutoCreateColumnsFromModel(true);
		this.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);

		// ((DefaultTableCellRenderer)this.getDefaultRenderer(String.class)).setFont(new Font("SansSerif", Font.BOLD, 10));

		JTableHeader header = getTableHeader();
		header.setDefaultRenderer(new HeaderRenderer(this));
		header.setFont(new Font("SansSerif", Font.BOLD, 10));

		setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

		/*
		// Figure out the row heights
		for (int i = 0; i < entries.size(); i++) {
			List<?> factors = (List)entries.get(i).get(HCAMetadata.FACTORS);
			int lines = factors.size();
			this.setRowHeight(i, lines*16);
		}
		*/

		// Set up selection listener
		getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent event) {
				int row = getSelectedRow();
				if (row != -1) {
					String acc = (String)getValueAt(row,0);
					hcaSource.getEntryFrame().enableButtons(acc, true);
				}
			}
		});

		// Set up double-click listener
		addMouseListener(new MouseAdapter() {
			public void mousePressed(MouseEvent mouseEvent) {
				JTable table = (JTable) mouseEvent.getSource();
				Point point = mouseEvent.getPoint();
				int row = table.rowAtPoint(point);
				if (mouseEvent.getClickCount() == 2 && table.getSelectedRow() != -1) {
					// Load the experiment
					String acc = (String)table.getValueAt(row,0);
					hcaSource.getEntryFrame().loadExperiment(acc);
					// TaskIterator tasks = new TaskIterator(new LoadExperimentTask(acc));
					// scNVManager.executeTasks(tasks);
				}
			}
		});
	}

	@Override
	public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
		Component returnComp = super.prepareRenderer(renderer, row, column);

		returnComp.setFont(new Font("SansSerif", Font.PLAIN, 10));
		if (!returnComp.getBackground().equals(getSelectionBackground())) {
			Color bg = (row % 2 == 0 ? alternateColor : Color.WHITE);
			returnComp.setBackground(bg);
			bg = null;
		}
		if (column == 2) {
			// Set the tooltip for this column
			int modelRow = convertRowIndexToModel(row);
			String ttText = createToolTip(getModel().getValueAt(modelRow, column).toString());
			((JLabel)returnComp).setToolTipText(ttText);
		} else {
			((JLabel)returnComp).setToolTipText(null);
		}
		return returnComp;
	}

	private String createToolTip(String tt) {
		return "<html><body><div style=\"width:300px\">"+tt+"</div></body></html>";
	}


}
