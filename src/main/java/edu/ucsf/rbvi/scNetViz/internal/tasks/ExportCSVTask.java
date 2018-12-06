package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVWriter;

// Tunable to choose experiment?

public class ExportCSVTask extends AbstractTask {
	final ScNVManager manager;
	final Matrix matrix;

	@Tunable(description="Output file: ",params="input=false")
	public File output = null;

	@Tunable(description="Delimiter: ")
	public ListSingleSelection<String> delimiter = new ListSingleSelection<String>(",","\t");

	public ExportCSVTask(final ScNVManager manager, Matrix matrix) {
		super();
		this.manager = manager;
		this.matrix = matrix;
		delimiter.setSelectedValue(",");
	}

	public void run(TaskMonitor monitor) {
		monitor.setTitle("Exporting CSV file to "+output.toString());
		CSVWriter.writeCSV(output, matrix, delimiter.getSelectedValue());
	}

}
