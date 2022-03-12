package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.List;

public interface Matrix {
	public String getMatrixType();
	public Class<?> getMatrixClass();
	public String toString();

	public int getNCols();
	public int getNRows();
	public int getNonZeroCount();
	public boolean isTransposed();
	public void setTranspose(boolean t);
  public int getHdrCols();
  public void setHdrCols(int hdrCols);
  public int getHdrRows();
  public void setHdrRows(int hdrRows);
  public int getRowKey();
  public void setRowKey(int keyRow);
  public int getColKey();
  public void setColKey(int keyCol);

	public List<String[]> getRowLabels();
	public List<String> getRowLabels(int col);
	public void setRowLabels(List<String[]> rLabels);
	public void setRowLabels(List<String> rLabels, int lbl);
	public void setRowLabels(String[] rLabels, int lbl);
	public void setRowLabel(String rLabel, int row, int col);
	public void setRowLabel(String[] rLabels, int row);

	public List<String[]> getColLabels();
	public List<String> getColLabels(int row);
	public void setColLabels(List<String[]> cLabels);
	public void setColLabels(String[] cLabels, int lbl);
	public void setColLabels(List<String> cLabels, int lbl);
	public void setColLabel(String cLabel, int row, int col);
	public void setColLabel(String[] cLabel, int col);

	public String getRowLabel(int row);
	public String getRowLabel(int row, int col);

	public String getColumnLabel(int col);
	public String getColumnLabel(int row, int col);

  public abstract Object getValue(int row, int col);
}


