package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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
	public static int LABEL_INDEX = 0;

	protected boolean transposed = false;

	protected int nRows;
	protected int nCols;
	protected int nonZeros;
  // The number of columns that contain row labels
  protected int hdrCols = 1;
  // The number of rows that contain column labels
  protected int hdrRows = 1; // NOTE: by default we don't include the first row

	protected List<String[]> rowLabels;
	protected List<String[]> colLabels;

  // The column header that contains the key information
	protected int columnKey = 0;
  // The row header that contains the key information
	protected int rowKey = LABEL_INDEX;

	protected String name;

	protected final ScNVManager manager;
	protected final CyTableFactory tableFactory;
	protected final CyTableManager tableManager;

	public SimpleMatrix(final ScNVManager manager) {
		this(manager, null, null);
	}

	public SimpleMatrix(final ScNVManager manager, 
	                    List<String[]> rowLabels, List<String[]> colLabels) {
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
  public int getHdrCols() { return hdrCols; }

  @Override
  public int getHdrRows() { return hdrRows; }

  @Override
  public void setHdrCols(int hdrCols) { this.hdrCols = hdrCols; }

  @Override
  public void setHdrRows(int hdrRows) { this.hdrRows = hdrRows; }

	@Override
	public int getNonZeroCount() { return nonZeros; }

	@Override
	public boolean isTransposed() { return transposed; }

	@Override
	public void setTranspose(boolean t) { transposed = t; }
  
  @Override
  public int getRowKey() { return rowKey; }
  
  @Override
  public void setRowKey(int rowKey) { this.rowKey = rowKey; }
  
  @Override
  public int getColKey() { return columnKey; }
  
  @Override
  public void setColKey(int colKey) { this.columnKey = colKey; }

	@Override
	public List<String[]> getRowLabels() { 
		return transposed ? colLabels : rowLabels;
	}

  @Override
  public List<String> getRowLabels(int col) {
    return getLabels(rowLabels, col);
  }

	@Override
	public void setRowLabels(List<String[]> rLabels) { rowLabels = rLabels; }

	@Override
	public void setRowLabel(String rLabel, int row, int col) { 
    if (rowLabels == null) {
      String[][] lbls = new String[nRows+hdrRows][];
      Arrays.fill(lbls, null);
      rowLabels = Arrays.asList(lbls);
    }

    if (rowLabels.size() <= row || rowLabels.get(row) == null)
      rowLabels.set(row, new String[hdrCols]);

    rowLabels.get(row)[col] = rLabel;
  }

	@Override
	public void setRowLabels(List<String> rLabels, int lbl) { 
    rowLabels = setLabels(rowLabels, rLabels.toArray(new String[0]), lbl);
  }

	@Override
	public void setRowLabels(String[] rLabels, int lbl) { 
    rowLabels = setLabels(rowLabels, rLabels, lbl);
  }

	@Override
	public void setRowLabel(String[] rLabels, int row) { 
    for (int col = 0; col < hdrCols; col++)
      setRowLabel(rLabels[col], row, col);
  }

	@Override
	public List<String[]> getColLabels() { 
		return transposed ? rowLabels : colLabels;
	}

  @Override
  public List<String> getColLabels(int row) {
    return getLabels(colLabels, row);
  }

	@Override
	public void setColLabels(List<String[]> cLabels) { colLabels = cLabels; }

	@Override
	public void setColLabel(String cLabel, int row, int col) { 
    if (colLabels == null) {
      String[][] lbls = new String[nCols+hdrCols][];
      Arrays.fill(lbls, null);
      colLabels = Arrays.asList(lbls);
    }

    if (colLabels.get(col) == null)
      colLabels.set(col, new String[hdrRows]);

    colLabels.get(col)[row] = cLabel;
  }

	@Override
	public void setColLabels(List<String> cLabels, int lbl) { 
    colLabels = setLabels(colLabels, cLabels.toArray(new String[0]), lbl);
  }

	@Override
	public void setColLabels(String[] cLabels, int lbl) { 
    colLabels = setLabels(colLabels, cLabels, lbl);
  }

	@Override
	public void setColLabel(String[] cLabels, int col) { 
    for (int row = 0; row < hdrRows; row++)
      setColLabel(cLabels[row], row, col);
  }

  private List<String[]> setLabels(List<String[]> lbls, String[] labels, int lbl) {
    if (lbls == null) {
      lbls = new ArrayList<>();
      for (String l: labels) {
        String[] lblArray = new String[hdrCols];
        lblArray[lbl] = l;
        lbls.add(lblArray);
      }
    } else {
      for (int i = 0; i < labels.length; i++) {
        lbls.get(i)[lbl] = labels[i];
      }
    }
    return lbls;
  }

	@Override
	public String getRowLabel(int row) {
    return getRowLabel(row,0);
  }

  @Override
  public String getRowLabel(int row, int col) {
		return transposed ? colLabels.get(row)[col] : rowLabels.get(row)[col];
	}

	@Override
	public String getColumnLabel(int col) {
    return getColumnLabel(col, 0);
  }

  @Override
  public String getColumnLabel(int col, int row) {
		String label = transposed ? rowLabels.get(col)[row]: colLabels.get(col)[row];
		return label;
	}

  public List<String> getLabels(List<String[]> labels, int index) {
    List<String>lbls = new ArrayList<>(labels.size());
    for (String[] lbl: labels)
      if (lbl != null)
        lbls.add(lbl[index]);
    return lbls;

  }
	public static List<String[]> getLabels(List<String[]> labelTable, int index, int hdrs, int extraHdrs, int hdr) {
    // System.out.println("getLabels: index="+index+", hdrs="+hdrs+", hdr="+hdr);
		List<String[]> labels = new ArrayList<>(labelTable.size()+extraHdrs);
    if (extraHdrs > 0) {
      for (int i = 0; i < extraHdrs; i++) {
        String[] lbl = new String[hdrs];
        lbl[hdr] = "";
			  labels.add(lbl);
      }
    }
		for (String[] row: labelTable) {
      String[] lbl = new String[hdrs];
      lbl[hdr] = row[index];
			labels.add(lbl);
		}
		return labels;
	}

	public static List<String[]> getLabels(List<String[]> labelTable, int hdrs, int hdr) {
		return getLabels(labelTable, LABEL_INDEX, hdrs, 0, hdr);
	}

}
