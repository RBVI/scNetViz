package edu.ucsf.rbvi.scNetViz.internal.view;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class DefaultExperimentTableModel extends SortableTableModel {
	final Experiment experiment;
	final DoubleMatrix matrix;

	public DefaultExperimentTableModel (final Experiment experiment) {
		super(1); // One header row by default
		this.experiment = experiment;
		if (experiment.getMatrix() instanceof DoubleMatrix)
			matrix = (DoubleMatrix)experiment.getMatrix();

		// Throw an exception?
		throw new RuntimeException("Experiment matrices must be of type double");
	}

	@Override
	public int getColumnCount() { return matrix.getNCols(); }

	@Override
	public int getSelectedRow() { return 0; }

	@Override
	public void setSelectedRow(int selectedRow) {  }

	@Override
	public String getColumnName(int column) {
		if (column == 0) {
			return matrix.isTransposed() ? "Barcodes" : "Genes";
		}
		if (columnIndex != null) {
			return matrix.getColumnLabel(columnIndex[column-1]);
		}
		return matrix.getColumnLabel(column-1);
	}

	@Override
	public int getRowCount() { 
		return matrix.getNRows(); 
	}

	@Override
	public Class getColumnClass(int column) {
		switch (column) {
			case 0:
				return String.class;
			default:
				return Double.class;
		}
	}

	@Override
	public Object getValueAt(int row, int column) {
		//System.out.println("getValueAt: "+row+","+column);
		switch (column) {
			case 0:
				return matrix.getRowLabel(row);
			default:
				double v;
			 	if (columnIndex != null)
					v	= matrix.getDoubleValue(row, columnIndex[column]);
				else
					v	= matrix.getDoubleValue(row, column);
				if (Double.isNaN(v)) return null;
				return new Double(v);
		}
	}

}
