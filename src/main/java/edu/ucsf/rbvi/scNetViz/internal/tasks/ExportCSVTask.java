package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;

import java.util.Arrays;

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
	Matrix matrix;

	enum DELIMITERS {
		COMMA(",",","),
		TAB("TAB", "\t");

		String text;
		String delim;

		DELIMITERS(String text, String delim) {
			this.text = text;
			this.delim = delim;
		}

		public String toString() { return text; }
		public String delim() { return delim; }
	};

	@Tunable(description="Output file: ",params="input=false")
	public File file = null;

	@Tunable(description="Delimiter: ")
	public ListSingleSelection<DELIMITERS> delimiter = 
		new ListSingleSelection<DELIMITERS>(Arrays.asList(DELIMITERS.values()));

	public ExportCSVTask(final ScNVManager manager) {
		super();
		this.manager = manager;
		delimiter.setSelectedValue(DELIMITERS.COMMA);
	}

	public ExportCSVTask(final ScNVManager manager, Matrix matrix) {
		super();
		this.manager = manager;
		this.matrix = matrix;
		delimiter.setSelectedValue(DELIMITERS.COMMA);
	}

	public void setMatrix(Matrix matrix) { this.matrix = matrix; }

	public void run(TaskMonitor monitor) {
		System.out.println("ExportCSVTask");
		monitor.setTitle("Exporting CSV file to "+file.toString());
		CSVWriter.writeCSV(file, matrix, delimiter.getSelectedValue().delim());
	}

}
