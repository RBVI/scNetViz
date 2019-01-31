package edu.ucsf.rbvi.scNetViz.internal.api;

public interface DoubleMatrix extends Matrix {
	public double[][] getDoubleMatrix(double missing);
	public double[][] getDoubleMatrix(double missing, boolean transpose);

	// Maybe this should be more general?
	public double[][] getDoubleMatrix(double missing, boolean transpose, boolean excludeControls);
	public double getDoubleValue(String rowLabel, String colLabel);
	public double getDoubleValue(int rowIndex, int colIndex);
	public default Class<?> getMatrixClass() { return Double.class; }
}
