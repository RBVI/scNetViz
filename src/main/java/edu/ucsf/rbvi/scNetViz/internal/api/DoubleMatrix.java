package edu.ucsf.rbvi.scNetViz.internal.api;

public interface DoubleMatrix extends Matrix {
	public double[][] getDoubleMatrix(double missing);
	public double getDoubleValue(String rowLabel, String colLabel);
	public double getDoubleValue(int rowIndex, int colIndex);
}
