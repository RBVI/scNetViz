package edu.ucsf.rbvi.scNetViz.internal.api;

public interface IntegerMatrix extends Matrix {
	public int[][] getIntegerMatrix(int missing);
	public int[][] getIntegerMatrix(int missing, boolean transpose);
	public int getIntegerValue(String rowLabel, String colLabel);
	public int getIntegerValue(int rowIndex, int colIndex);
	public default Class<?> getMatrixClass() { return Integer.class; }
}
