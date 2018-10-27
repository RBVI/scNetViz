package edu.ucsf.rbvi.scNetViz.internal.api;

import java.util.List;

public interface Matrix {
	public String getMatrixType();
	public String toString();

	public int getNCols();
	public int getNRows();
	public int getNonZeroCount();
	public boolean isTransposed();
	public void setTranspose(boolean t);

	public List<String> getRowLabels();
	public void setRowLabels(List<String> rLabels);

	public List<String> getColLabels();
	public void setColLabels(List<String> cLabels);

	public String getRowLabel(int row);

	public String getColumnLabel(int col);
}


