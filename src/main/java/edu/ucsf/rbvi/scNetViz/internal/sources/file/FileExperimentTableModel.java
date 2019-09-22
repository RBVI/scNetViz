package edu.ucsf.rbvi.scNetViz.internal.sources.file;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public class FileExperimentTableModel extends SortableTableModel {
	final ScNVManager manager;
	final FileExperiment fileExperiment;
	final MatrixMarket matrixMarket;

	public FileExperimentTableModel (final ScNVManager manager, 
	                                final FileExperiment experiment) {
		super(1);
		this.manager = manager;
		this.fileExperiment = experiment;
		this.matrixMarket = (MatrixMarket)experiment.getMatrix();
	}

	@Override
	public int getColumnCount() { return matrixMarket.getNCols()+hdrCols; }

	@Override
	public int getSelectedRow() { return 0; }

	@Override
	public void setSelectedRow(int selectedRow) {  }

	@Override
	public String getColumnName(int column) {
		if (column == 0) {
			return matrixMarket.isTransposed() ? "Barcodes" : "Genes";
		}
		if (columnIndex != null) {
			return matrixMarket.getColumnLabel(columnIndex[column-1]);
		}
		return matrixMarket.getColumnLabel(column);
	}

	@Override
	public int getRowCount() { 
		return matrixMarket.getNRows(); 
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
		// System.out.println("getValueAt: "+row+","+column);
		switch (column) {
			case 0:
				return matrixMarket.getRowLabel(row);
			default:
				double v;
			 	if (columnIndex != null)
					v	= matrixMarket.getDoubleValue(row, columnIndex[column-1]);
				else
					v = matrixMarket.getDoubleValue(row, column-1);
				if (Double.isNaN(v)) return null;
				return new Double(v);
		}
	}

}
