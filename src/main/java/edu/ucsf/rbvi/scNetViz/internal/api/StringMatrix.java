package edu.ucsf.rbvi.scNetViz.internal.api;

public interface StringMatrix extends Matrix {
	public String[][] getStringMatrix();
	public String[][] getStringMatrix(boolean tranpose);
	public String[][] getStringMatrix(boolean tranpose, boolean excludeControls);
	public String getValue(String rowLabel, String colLabel);
	public String getValue(int rowIndex, int colIndex);
	public default Class<?> getMatrixClass() { return String.class; }
}
