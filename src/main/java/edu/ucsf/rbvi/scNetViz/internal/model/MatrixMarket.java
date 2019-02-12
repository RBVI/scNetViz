package edu.ucsf.rbvi.scNetViz.internal.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;
import org.cytoscape.model.CyTableFactory;
import org.cytoscape.model.CyTableManager;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.MatrixUtils;

public class MatrixMarket extends SimpleMatrix implements DoubleMatrix, IntegerMatrix {
	public static String HEADER = "%%MatrixMarket";
	public static String COMMENT = "%";
	public static String delimiter = null;
	public static int LABEL_INDEX = 1;

	public static enum MTXOBJECT {
		MATRIX("matrix"),
		// DGRAPH("directed graph"),
		VECTOR("vector");
	
		String strOType;
		MTXOBJECT(String type) {
			this.strOType = type;
		}

		public String toString() { return strOType; }

		public static MTXOBJECT getEnum(String str) {
			for (MTXOBJECT obj: MTXOBJECT.values()) {
				if (str.toLowerCase().equals(obj.toString()))
					return obj;
			}
			return null;
		}
	}


	public enum MTXFORMAT {
		COORDINATE("coordinate"),
		ARRAY("array");
	
		String strFormat;
		MTXFORMAT(String format) {
			this.strFormat = format;
		}

		public String toString() { return strFormat; }

		public static MTXFORMAT getEnum(String str) {
			for (MTXFORMAT obj: MTXFORMAT.values()) {
				if (str.toLowerCase().equals(obj.toString()))
					return obj;
			}
			return null;
		}
	}

	public enum MTXTYPE {
		REAL("real"),
		COMPLEX("complex"),
		INTEGER("integer"),
		PATTERN("pattern");
	
		String strType;
		MTXTYPE(String type) {
			this.strType = type;
		}

		public String toString() { return strType; }

		public static MTXTYPE getEnum(String str) {
			for (MTXTYPE obj: MTXTYPE.values()) {
				if (str.toLowerCase().equals(obj.toString()))
					return obj;
			}
			return null;
		}
	}

	public enum MTXSYMMETRY {
		GENERAL("general"),
		SYMTXETRIC("symmetric"),
		SKEW("skew-symmetric"),
		HERMITIAN("hermitian");

		String strSymmetry;
		MTXSYMMETRY(String sym) {
			this.strSymmetry = sym;
		}

		public String toString() { return strSymmetry; }

		public static MTXSYMMETRY getEnum(String str) {
			for (MTXSYMMETRY obj: MTXSYMMETRY.values()) {
				if (str.toLowerCase().equals(obj.toString()))
					return obj;
			}
			return null;
		}
	}


	// Information about the matrix
	MTXOBJECT objectType;
	MTXFORMAT format;
	MTXTYPE type;
	MTXSYMMETRY sym;

	private List<String> comments;
	private Map<String, Integer> rowMap;
	private Map<String, Integer> colMap;

	// We only support real and integer at this point
	private int[][] intMatrix;
	private double[][] doubleMatrix;

	private int[] colIndex;

	private List<String[]> rowTable;
	private List<String[]> colTable;
	private BitSet controls;

	public MatrixMarket(final ScNVManager manager) {
		this(manager, null, null);
	}

	public MatrixMarket(final ScNVManager manager, 
	                    List<String> rowLabels, List<String> colLabels) {
		super(manager, rowLabels, colLabels);
		rowMap = new HashMap<>();
		colMap = new HashMap<>();
	}

	public MTXFORMAT getFormat() {
		return format;
	}

	public MTXOBJECT getObjectType() {
		return objectType;
	}

	public MTXTYPE getType() {
		return type;
	}

	public MTXSYMMETRY getSymmetry() {
		return sym;
	}

	public BitSet findControls() {
		if (controls == null)
			controls = MatrixUtils.findControls(rowLabels, MatrixUtils.CONTROL_PREFIX);
		return controls;
	}

	public boolean isControl(int row) {
		if (controls == null)
			controls = MatrixUtils.findControls(rowLabels, MatrixUtils.CONTROL_PREFIX);

		return controls.get(row);
	}

	@Override
	public String getMatrixType() { return "MatrixMarket"; }

	@Override
	public Class<?> getMatrixClass() {
		if (type == MTXTYPE.INTEGER)
			return Integer.class;
		return Double.class;
	}

	public void setRowTable(List<String[]> rowTable) {
		this.rowTable = rowTable;
		if (rowTable != null)
			rowLabels = getLabels(rowTable);
	}

	public void setColumnTable(List<String[]> colTable) {
		this.colTable = colTable;
		if (colTable != null) {
			// System.out.println("Column tables has "+colTable.size()+" entries");
			colLabels = getLabels(colTable);
		}
	}

	public void readMTX(TaskMonitor taskMonitor, File mmInputName) throws FileNotFoundException, IOException {
		FileInputStream inputStream = new FileInputStream(mmInputName);
		readMTX(taskMonitor, inputStream, mmInputName.getName());
	}

	public void readMTX(TaskMonitor taskMonitor, InputStream stream, String mmInputName) throws FileNotFoundException, IOException {
		this.name = mmInputName;

		//System.out.println("Reading "+name);

		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		// Read the first line
		String header = reader.readLine();
		parseHeader(header);

		//System.out.println("Got header");

		// Now, read until we find the dimensions
		comments = new ArrayList<>();
		String line = reader.readLine();
		while (line.startsWith(COMMENT) || line.trim().length() == 0) {
			if (line.startsWith(COMMENT))
				comments.add(line);
			line = reader.readLine();
		}

		// at this point, line should have our dimensions.
		// if we have a coordinate, we want three values, otherwise we want two
		String[] dims = line.split("\\s+");
		nRows = Integer.parseInt(dims[0]);
		nCols = Integer.parseInt(dims[1]);
		//System.out.println("nRows = "+nRows);
		//System.out.println("nCols = "+nCols);
		if (format == MTXFORMAT.ARRAY) {
			if (type == MTXTYPE.REAL) {
				doubleMatrix = new double[nRows][nCols];
			} else if (type == MTXTYPE.INTEGER) {
				intMatrix = new int[nRows][nCols];
			}

			for (int col = 0; col < nCols; col++) {
				for (int row = 0; row < nRows; row++) {
					readArrayLine(row, col, reader);
				}
			}
		} else if (format == MTXFORMAT.COORDINATE) {
			nonZeros = Integer.parseInt(dims[2]);
			// System.out.println("nonZeros = "+nonZeros);
			colIndex = new int[nCols+1];
			Arrays.fill(colIndex, -1);
			colIndex[nCols] = nonZeros; // Point to the end
			if (type == MTXTYPE.INTEGER) {
				intMatrix = new int[nonZeros][3];
			} else if (type == MTXTYPE.REAL) {
				intMatrix = new int[nonZeros][2]; // indices
				doubleMatrix = new double[nonZeros][1]; // actual data
			}
			for (int index = 0; index < nonZeros; index++) {
				readCoordinateLine(index, reader);
			}
		}
	}

	public void mergeTable(CyTable table, String mergeColumn) {
	}

	public CyTable makeTable(String name, boolean addTable) {
		// Create the table
		CyTable table;
		if (transposed) {
			table = tableFactory.createTable(name, "Cells", String.class, true, false);
			createColumns(table, rowLabels);
			// System.out.println("Created "+rowLabels.size()+" columns");
			createRows(table, colLabels, rowLabels);
			// System.out.println("Created "+colLabels.size()+" rows");
		} else {
			table = tableFactory.createTable(name, "Genes", String.class, true, false);
			createColumns(table, colLabels);
			// System.out.println("Created columns");
			createRows(table, rowLabels, colLabels);
			// System.out.println("Created rows");
		}

		if (addTable) {
			// System.out.println("Adding table to table manager");
			tableManager.addTable(table);
			// System.out.println("...done...");
		}
		return table;
	}

	public int[][] getIntegerMatrix(int missing) {
		return getIntegerMatrix(missing, false, false);
	}

	public int[][] getIntegerMatrix(int missing, boolean transpose) {
		return getIntegerMatrix(missing, transpose, false);
	}

	public int[][] getIntegerMatrix(int missing, boolean transpose, boolean excludeControls) {
		BitSet excludeRows = null;
		if (excludeControls) {
			excludeRows = findControls();
		}

		// Are we tranposing or not?
		if (transpose && transposed)
			transpose = false;
		else
			transpose = transpose | transposed;

		if (format == MTXFORMAT.ARRAY) {
			if (type == MTXTYPE.INTEGER)
				if (!transpose && !excludeControls)
					return intMatrix;
				else {
					int[][] newArray = getIntegerMatrix(transpose, excludeRows);

					int newRow = 0;
					for (int row = 0; row < nRows; row++) {
						if (excludeRows != null && excludeRows.get(row))
							continue;
						for (int col = 0; col < nCols; col++) {
							if (transposed)
								newArray[col][newRow++] = intMatrix[col][row];
							else
								newArray[newRow++][col] = intMatrix[row][col];
						}
					}
					return newArray;
				}

			if (type == MTXTYPE.REAL) {
				int[][] newArray = getIntegerMatrix(transpose, excludeRows);

				int newRow = 0;
				for (int row = 0; row < nRows; row++) {
					if (excludeRows != null && excludeRows.get(row))
						continue;
					for (int col = 0; col < nCols; col++) {
						if (transposed)
							newArray[col][newRow++] = (int)Math.round(doubleMatrix[col][row]);
						else
							newArray[newRow++][col] = (int)Math.round(doubleMatrix[row][col]);
					}
				}
				return newArray;
			}
		} else if (format == MTXFORMAT.COORDINATE) {
			int[][] newArray = getIntegerMatrix(transpose, excludeRows);

			int maxFill = getNRows();
			if (excludeRows != null) maxFill = maxFill-excludeRows.cardinality();
			for (int row = 0; row < maxFill; row++) {
				Arrays.fill(newArray[row], missing);
			}
			int newRow = 0;
			for (int index = 0; index < nonZeros; index++) {
				int row = intMatrix[index][0];
				if (excludeRows != null && excludeRows.get(row))
					continue;
				int col = intMatrix[index][1];
				row = newRow++;
				if (transposed) { int rtmp = row; row = col; col = rtmp; }
				if (type == MTXTYPE.INTEGER)
					newArray[row][col] = intMatrix[index][2];
				else if (type == MTXTYPE.REAL)
					newArray[row][col] = (int)Math.round(doubleMatrix[index][0]);
			}
			return newArray;
		}
		return null;
	}

	public int getIntegerValue(int row, int col) {
		if (transposed) { int rtmp = row; row = col; col = rtmp; }
		if (format == MTXFORMAT.ARRAY) {
			if (type == MTXTYPE.REAL)
				return (int)Math.round(doubleMatrix[row][col]);
			else if (type == MTXTYPE.INTEGER)
				return intMatrix[row][col];
		} else {
			// MTXFORMAT.COORDINATE
			int index = findIndex(row, col);
			if (index >= 0) {
				if (type == MTXTYPE.REAL)
					return (int)Math.round(doubleMatrix[index][0]);
				else if (type == MTXTYPE.INTEGER)
					return intMatrix[index][2];
			}
		}
		return Integer.MIN_VALUE;
	}

	public double getDoubleValue(int row, int col) {
		if (transposed) { int rtmp = row; row = col; col = rtmp; }
		if (format == MTXFORMAT.ARRAY) {
			if (type == MTXTYPE.REAL)
				return doubleMatrix[row][col];
			else if (type == MTXTYPE.INTEGER)
				return (double)intMatrix[row][col];
		} else {
			// MTXFORMAT.COORDINATE
			int index = findIndex(row, col);
			if (index >= 0) {
				if (type == MTXTYPE.REAL) {
					return doubleMatrix[index][0];
				} else if (type == MTXTYPE.INTEGER) {
					return (double)intMatrix[index][2];
				}
			}
		}
		return Double.NaN;
	}

	public double[][] getDoubleMatrix(double missing) {
		return getDoubleMatrix(missing, false);
	}

	public double[][] getDoubleMatrix(double missing, boolean transpose) {
		return getDoubleMatrix(missing, transpose, false);
	}

	public double[][] getDoubleMatrix(double missing, boolean transpose, boolean excludeControls) {
		BitSet excludeRows = null;
		if (excludeControls) {
			excludeRows = findControls();
		}
		// Are we tranposing or not?
		if (transpose && transposed)
			transpose = false;
		else
			transpose = transpose | transposed;
		//
		if (format == MTXFORMAT.ARRAY) {
			if (type == MTXTYPE.REAL) {
				if (!transpose && !excludeControls) {
					return doubleMatrix;
				} else {
					double[][] newArray = getDoubleMatrix(transpose, excludeRows);

					int newRow = 0;
					for (int row = 0; row < nRows; row++) {
						if (excludeRows != null && excludeRows.get(row))
							continue;
						for (int col = 0; col < nCols; col++) {
							if (transpose)
								newArray[col][newRow++] = (double)intMatrix[col][row];
							else
								newArray[newRow++][col] = (double)intMatrix[col][row];
						}
					}
					return newArray;
				}
			}

			if (type == MTXTYPE.INTEGER) {
				double[][] newArray = getDoubleMatrix(transpose, excludeRows);

				int newRow = 0;
				for (int row = 0; row < nRows; row++) {
					if (excludeRows != null && excludeRows.get(row))
						continue;
					for (int col = 0; col < nCols; col++) {
						if (transpose)
							newArray[col][newRow++] = (double)intMatrix[col][row];
						else
							newArray[newRow++][col] = (double)intMatrix[col][row];
					}
				}
				return newArray;
			}
		} else if (format == MTXFORMAT.COORDINATE) {
			double[][] newArray = getDoubleMatrix(transpose, excludeRows);

			int maxFill = getNRows();
			for (int row = 0; row < newArray.length; row++) {
				Arrays.fill(newArray[row], missing);
			}
			int skippedRows = 0;
			int lastRow = 0;
			for (int index = 0; index < nonZeros; index++) {
				int row = intMatrix[index][0];
				if (row < lastRow) {
					// We've wrapped
					skippedRows = 0;
				}
				lastRow = row;
				if (excludeRows != null && excludeRows.get(row)) {
					skippedRows++;
					continue;
				} else if (excludeRows != null) {
					row = row - skippedRows;
				}
				int col = intMatrix[index][1];
				// System.out.println("Row = "+row+", Col = "+col);
				if (transpose) { int rtmp = row; row = col; col = rtmp; }
				if (type == MTXTYPE.INTEGER)
					newArray[row][col] = (double)intMatrix[index][2];
				else if (type == MTXTYPE.REAL) {
					// System.out.println("newArray["+row+"]["+col+"] = "+doubleMatrix[index][0]+", skippedRows = "+skippedRows);
					newArray[row][col] = doubleMatrix[index][0];
				}
			}
			return newArray;
		}
		return null;
	}

	public int getIntegerValue(String rowLabel, String colLabel) {
		int row = rowLabels.indexOf(rowLabel);
		int col = colLabels.indexOf(colLabel);
		if (transposed) {
			int tmp = row; row = col; col = tmp;
		}
		int v = getIntegerValue(row+1, col+1);
		// System.out.println("Value for "+rowLabel+":"+row+","+colLabel+":"+col+" = "+v);
		return v;
	}

	public double getDoubleValue(String rowLabel, String colLabel) {
		int row = rowLabels.indexOf(rowLabel);
		int col = colLabels.indexOf(colLabel);
		if (transposed) {
			int tmp = row; row = col; col = tmp;
		}
		double v = getDoubleValue(row+1, col+1);
		// System.out.println("Value for "+rowLabel+":"+row+","+colLabel+":"+col+" = "+v);
		return v;
	}

	private void parseHeader(String header) throws IOException {
		if (!header.startsWith(HEADER))
			throw new IOException("File doesn't start with appropriate header");

		String[] headerArgs = header.split("\\s+");

		objectType = MTXOBJECT.getEnum(headerArgs[1]);
		if (objectType == null)
			throw new IOException("Illegal or unsupported object type: "+headerArgs[1]);
		if (objectType == MTXOBJECT.VECTOR)
			throw new IOException("Vector objects are not supported at this time");

		format = MTXFORMAT.getEnum(headerArgs[2]);
		if (format == null)
			throw new IOException("Illegal or unsupported format: "+headerArgs[2]);

		type = MTXTYPE.getEnum(headerArgs[3]);
		if (type == null)
			throw new IOException("Illegal or unsupported type: "+headerArgs[3]);
		if (type == MTXTYPE.COMPLEX || type == MTXTYPE.PATTERN)
			throw new IOException("Complex and pattern types are not supported at this time");

		sym = MTXSYMMETRY.getEnum(headerArgs[4]);
		if (sym == null)
			throw new IOException("Illegal or unsupported symmetry: "+headerArgs[4]);
		if (sym == MTXSYMMETRY.HERMITIAN)
			throw new IOException("Hermitian symmetry is not supported at this time");
	}

	private void readArrayLine(int row, int col, BufferedReader reader) throws IOException {
		String line = reader.readLine();
		// Skip over blank lines
		while (line.trim().length() == 0)
			line = reader.readLine();
		if (type == MTXTYPE.INTEGER) {
			intMatrix[row][col] = Integer.parseInt(line.trim());
		} else if (type == MTXTYPE.REAL) {
			doubleMatrix[row][col] = Double.parseDouble(line.trim());
		}
	}

	private void readCoordinateLine(int index, BufferedReader reader) throws IOException {
		String line = reader.readLine();
		// Skip over blank lines
		while (line.trim().length() == 0)
			line = reader.readLine();

		// We should have exactly three value: row, col, value
		String[] vals = line.split("\\s+");

		// Subtract 1 to make everything 0 relative
		intMatrix[index][0] = Integer.parseInt(vals[0])-1;
		intMatrix[index][1] = Integer.parseInt(vals[1])-1;
		if (colIndex[intMatrix[index][1]] < 0) {
			colIndex[intMatrix[index][1]] = index;
		}
		if (type == MTXTYPE.INTEGER) {
			intMatrix[index][2] = Integer.parseInt(vals[2]);
		} else if (type == MTXTYPE.REAL) {
			doubleMatrix[index][0] = Double.parseDouble(vals[2]);
			// System.out.println("value["+intMatrix[index][0]+"]["+intMatrix[index][1]+"] = "+doubleMatrix[index][0]);
		}
	}

	private void createColumns(CyTable table, List<String> labels) {
		for (String lbl: labels) {
			String columnLabel = "MTX::"+lbl;
			table.createColumn(columnLabel, Double.class, true);
		}
	}

	private void createRows(CyTable table, List<String> rowLabels, List<String> colLabels) {
		if (format == MTXFORMAT.COORDINATE) {
			for (int index = 0; index < nonZeros; index++) {
				int row = intMatrix[index][0];
				int col = intMatrix[index][1];
				if (transposed) { int rtmp = row; row = col; col = rtmp; }
				String rowLabel = rowLabels.get(row);
				CyRow cyRow = table.getRow(rowLabel);
				String colLabel = "MTX::"+colLabels.get(col);
				if (type == MTXTYPE.REAL) {
					double v = doubleMatrix[index][0];
					cyRow.set(colLabel, Double.valueOf(v));
				} else if (type == MTXTYPE.INTEGER) {
					int v = intMatrix[index][2];
					cyRow.set(colLabel, Double.valueOf(v));
				}
			}
		} else if (format == MTXFORMAT.ARRAY) {
			for (int row = 0; row < rowLabels.size(); row++) {
				String rowLabel = rowLabels.get(row);
				CyRow cyRow = table.getRow(rowLabel);
				for (int col = 0; col < colLabels.size(); col++) {
					String colLabel = "MTX::"+colLabels.get(col);
					int r = row;
					int c = col;
					if (transposed) { int rtmp = r; r = c; c = rtmp; }
					if (type == MTXTYPE.REAL) {
						cyRow.set(colLabel, Double.valueOf(doubleMatrix[r][c]));
					} else if (type == MTXTYPE.INTEGER) {
						cyRow.set(colLabel, Double.valueOf(intMatrix[r][c]));
					}
				}
			}
		}
	}

	private int findIndex(int row, int col) {
		// System.out.println("findIndex("+row+","+col+")");
		int index1 = colIndex[col];
		// System.out.println("index1="+index1);
		if (index1 < 0) return Integer.MIN_VALUE;

		int index2 = colIndex[++col]-1;
		// System.out.println("index2="+index2);

		int indexMid = index1+(index2-index1)/2;
		int index = -1;
		while (index < 0) {
			index = compare(index1, index2, indexMid, row);
			// System.out.println("compare = "+index);
			if (index == -3) {
				return -1;
			} else if (index == -2) {
				index2 = indexMid;
				indexMid = index1+(index2-index1)/2;
				if (index2 == indexMid)
					return -1;
			} else if (index == -1) {
				index1 = indexMid;
				indexMid = index1+(index2-index1)/2;
				if (index1 == indexMid)
					return -1;
			}
		}
		//System.out.println("index="+index);
		return index;
	}

	private int compare(int index1, int index2, int indexMid, int row) {
		/*
		System.out.println("index1="+index1+", index2="+index2+", indexMid="+indexMid+", row = "+row);
		System.out.println("intMatrix[index1][0] = "+intMatrix[index1][0]);
		System.out.println("intMatrix[indexMid][0] = "+intMatrix[indexMid][0]);
		System.out.println("intMatrix[index2][0] = "+intMatrix[index2][0]);
		*/
		if (intMatrix[index1][0] > row || intMatrix[index2][0] < row)
			return -3;
		if (intMatrix[index1][0] < row && intMatrix[indexMid][0] > row)
			return -2;
		if (intMatrix[indexMid][0] < row && intMatrix[index2][0] > row)
			return -1;
		if (intMatrix[index1][0] == row) return index1;
		if (intMatrix[indexMid][0] == row) return indexMid;
		return index2;
	}

	private double[][] getDoubleMatrix(boolean transpose, BitSet excludeRows) {
		int excludeSize = 0;
		if (excludeRows != null)
			excludeSize = excludeRows.cardinality();
		if (transpose)
			return new double[nCols][nRows-excludeSize];
		else
			return new double[nRows-excludeSize][nCols];
	}

	private int[][] getIntegerMatrix(boolean transpose, BitSet excludeRows) {
		int excludeSize = 0;
		if (excludeRows != null)
			excludeSize = excludeRows.cardinality();
		if (transpose)
			return new int[nCols][nRows-excludeSize];
		else
			return new int[nRows-excludeSize][nCols];
	}
}
