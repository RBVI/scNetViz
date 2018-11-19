package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.view;

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
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;

import org.cytoscape.application.swing.CySwingApplication;
import org.cytoscape.application.swing.CytoPanel;
import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.CyTable;
import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAExperiment;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ProcessAllTask;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ShowExperimentTableTask;
import edu.ucsf.rbvi.scNetViz.internal.view.HeaderRenderer;

public class GXAEntryTable extends JTable {
	final ScNVManager scNVManager;
	final GXASource gxaSource;
	List<Metadata> entries;

	static String[] columnNames = {"Loaded", "Experiment", "Assays", "Comparisons", "Organisms", "Experimental Variables"};

	static Color alternateColor = new Color(234,255,234);

	public GXAEntryTable (final ScNVManager scNVManager, final GXASource gxaSource, final GXAEntryTableModel tableModel) {
		super(tableModel);
		this.scNVManager = scNVManager;
		this.gxaSource = gxaSource;

		entries = gxaSource.getMetadata();

		this.setAutoCreateRowSorter(true);
		this.setAutoCreateColumnsFromModel(true);
		this.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);

		// ((DefaultTableCellRenderer)this.getDefaultRenderer(String.class)).setFont(new Font("SansSerif", Font.BOLD, 10));

		JTableHeader header = getTableHeader();
		header.setDefaultRenderer(new HeaderRenderer(this));
		header.setFont(new Font("SansSerif", Font.BOLD, 10));

		// Figure out the row heights
		for (int i = 0; i < entries.size(); i++) {
			List<?> factors = (List)entries.get(i).get(GXAMetadata.FACTORS);
			int lines = factors.size();
			this.setRowHeight(i, lines*16);
		}

		addMouseListener(new MouseAdapter() {
			public void mousePressed(MouseEvent mouseEvent) {
				JTable table = (JTable) mouseEvent.getSource();
				Point point = mouseEvent.getPoint();
				int row = table.rowAtPoint(point);
				if (mouseEvent.getClickCount() == 2 && table.getSelectedRow() != -1) {
					// Load the experiment
					String acc = (String)table.getValueAt(row,0);
					TaskIterator tasks = new TaskIterator(new LoadExperimentTask(acc));
					scNVManager.executeTasks(tasks);
				}
			}
		});
	}

	@Override
	public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
		Component returnComp = super.prepareRenderer(renderer, row, column);

		if (!returnComp.getBackground().equals(getSelectionBackground())) {
			Color bg = (row % 2 == 0 ? alternateColor : Color.WHITE);
			returnComp.setBackground(bg);
			returnComp.setFont(new Font("SansSerif", Font.PLAIN, 10));
			bg = null;
		}
		return returnComp;
	}

	public class LoadExperimentTask extends AbstractTask {
		String accession;

		LoadExperimentTask(String acc) {
			this.accession = acc;
		}

		@Override
		public void run(TaskMonitor taskMonitor) {
			taskMonitor.setTitle("Loading Single Cell Expression Atlas experiment "+accession);
			GXAExperiment experiment = (GXAExperiment)gxaSource.getExperiment(accession, taskMonitor);
			if (experiment != null)
				scNVManager.addExperiment(accession,experiment);
			gxaSource.showEntriesTable(false);
			experiment.fetchClusters(taskMonitor);
			experiment.fetchDesign(taskMonitor);

			// Now, depending on the setting of our checkbox, either show
			// the experiments table or process the data
			if (gxaSource.getEntryFrame().isLoadOnly()) {
				// Create the Experiment table and show it
				insertTasksAfterCurrentTask(new ShowExperimentTableTask(scNVManager, experiment));
			} else {
				// Create the autoprocess task and append it
				insertTasksAfterCurrentTask(new ProcessAllTask(scNVManager, experiment));
			}
		}

		@ProvidesTitle
		public String getTitle() {return "Load Single Cell Expression Atlas Experiment "+accession;}
	}

}
