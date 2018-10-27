package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.CyTableFactory;
import org.cytoscape.model.CyTableManager;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;

public abstract class SimpleMatrix implements Matrix {
	public static int LABEL_INDEX = 1;

	protected boolean transposed = false;

	protected int nRows;
	protected int nCols;
	protected int nonZeros;

	protected List<String> rowLabels;
	protected List<String> colLabels;

	protected String name;

	protected final ScNVManager manager;
	protected final CyTableFactory tableFactory;
	protected final CyTableManager tableManager;

	public SimpleMatrix(final ScNVManager manager) {
		this(manager, null, null);
	}

	public SimpleMatrix(final ScNVManager manager, 
	                    List<String> rowLabels, List<String> colLabels) {
		this.manager = manager;
		this.rowLabels = rowLabels;
		this.colLabels = colLabels;
		this.tableFactory = manager.getService(CyTableFactory.class);
		this.tableManager = manager.getService(CyTableManager.class);
	}

	@Override
	public String toString() { return name; }

	@Override
	public int getNCols() { return transposed ? nRows : nCols; }

	@Override
	public int getNRows() { return transposed ? nCols : nRows; }

	@Override
	public int getNonZeroCount() { return nonZeros; }

	@Override
	public boolean isTransposed() { return transposed; }

	@Override
	public void setTranspose(boolean t) { transposed = t; }

	@Override
	public List<String> getRowLabels() { 
		return transposed ? colLabels : rowLabels;
	}

	@Override
	public void setRowLabels(List<String> rLabels) { rowLabels = rLabels; }

	@Override
	public List<String> getColLabels() { 
		return transposed ? rowLabels : colLabels;
	}

	@Override
	public void setColLabels(List<String> cLabels) { colLabels = cLabels; }

	@Override
	public String getRowLabel(int row) {
		return transposed ? colLabels.get(row) : rowLabels.get(row);
	}

	@Override
	public String getColumnLabel(int col) {
		return transposed ? rowLabels.get(col) : colLabels.get(col);
	}

	public static List<String> getLabels(List<String[]> labelTable) {
		List<String> labels = new ArrayList<>(labelTable.size());
		for (String[] row: labelTable) {
			labels.add(row[LABEL_INDEX]);
		}
		return labels;
	}

}
