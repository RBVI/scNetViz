package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.StringMatrix;

public class CSVWriter {

	public static void writeCSV(OutputStream output, List<String[]> table, int nCols) throws IOException {
		for (String[] line: table) {
			if (nCols < 0) nCols = line.length;
			int count = nCols;
			for (String txt: line) {
				output.write(txt.getBytes());
				if (count-- > 1)
					output.write("\t".getBytes());
				else
					break;
			}
			output.write("\n".getBytes());
		}
	}

	public static void writeCSV(OutputStream output, List<String[]> table) throws IOException {
		writeCSV(output, table, -1);
		/*
		for (String[] line: table) {
			int count = line.length;
			for (String txt: line) {
				output.write(txt.getBytes());
				if (count-- > 1)
					output.write("\t".getBytes());
			}
			output.write("\n".getBytes());
		}
		*/
	}

	public static void writeCSV(File file, List<String[]> table) throws IOException {
		BufferedWriter output = new BufferedWriter(new FileWriter(file));
		for (String[] line: table) {
			int count = line.length;
			for (String txt: line) {
				output.write(txt);
				if (count-- > 1)
					output.write("\t");
			}
			output.write("\n");
		}
		output.close();
	}

	public static void writeCSV(File file, Matrix matrix, String delimiter) {
		BufferedWriter output = null;
		try {
			output = new BufferedWriter(new FileWriter(file));
			List<String[]> rowLabels = matrix.getRowLabels();
			// Get the column labels
			List<String[]> colLabels = matrix.getColLabels();
      int hdrRows = matrix.getHdrRows();

      for (int hdrRow = 0; hdrRow < hdrRows; hdrRow++) {
			  for (int i = 0; i < colLabels.size()-1; i++) {
				  output.write(quote(colLabels.get(i)[hdrRow], delimiter)+delimiter);
        }
        if (colLabels.get(colLabels.size()-1) == null)
          output.write("\n");
        else
          output.write(quote(colLabels.get(colLabels.size()-1)[hdrRow], delimiter)+"\n");
			}

			// Get the matrix
			if (matrix.getMatrixClass() == Double.class) {
				for (int row = 0; row < rowLabels.size(); row++) {
					writeDoubleRow(output, rowLabels.get(row), matrix, row, delimiter);
				}
			} else if (matrix.getMatrixClass() == Integer.class) {
				for (int row = 0; row < rowLabels.size(); row++) {
					writeIntegerRow(output, rowLabels.get(row), matrix, row, delimiter);
				}
			} else if (matrix.getMatrixClass() == String.class) {
				for (int row = 0; row < rowLabels.size(); row++) {
					writeStringRow(output, rowLabels.get(row), matrix, row, delimiter);
				}
      }

			output.close();
		} catch (Exception fnf) {
			fnf.printStackTrace();
		}
	}

	private static String quote(String str, String delimiter) {
		if (delimiter.equals("\t"))
			return str;
		return "\""+str+"\"";
	}

	private static void writeStringRow(BufferedWriter output, String[] rowLabel, 
	                                   Matrix mat, int row, String delimiter) throws IOException {
		int nCols = mat.getNCols();
    System.out.println("writeStringRow: nCols = "+nCols);
    System.out.println("writeStringRow: rowsLabel.length = "+rowLabel.length);
		for (int hdr=0; hdr < rowLabel.length; hdr++) {
		  output.write(quote(rowLabel[hdr], delimiter)+delimiter);
    }
		for (int col = 0; col < nCols-1; col++) {
      String value = (String)mat.getValue(row, col);
      if (value == null)
				output.write(delimiter);
      else
        output.write(quote(value, delimiter)+delimiter);
		}
    String value = (String)mat.getValue(row, nCols-1);
    if (value == null)
      output.write("\n");
    else
      output.write(quote(value, delimiter)+"\n");
	}

	private static void writeIntegerRow(BufferedWriter output, String[] rowLabel, 
	                                    Matrix mat, int row, String delimiter) throws IOException {
		int nCols = mat.getNCols();
    System.out.println("writeIntegerRow: nCols = "+nCols);
    System.out.println("writeIntegerRow: rowsLabel.length = "+rowLabel.length);
		for (int hdr=0; hdr < rowLabel.length; hdr++) {
		  output.write(quote(rowLabel[hdr], delimiter)+delimiter);
    }
		for (int col = 0; col < nCols-1; col++) {
			Integer value = (Integer)mat.getValue(row, col);
			if (value == null || value == Integer.MIN_VALUE)
				output.write(delimiter);
			else
				output.write(value+delimiter);
		}
		Integer value = (Integer)mat.getValue(row, nCols-1);
    // System.out.println("Value["+row+","+(nCols-1)+"] = "+value);
		if (value == null || value == Integer.MIN_VALUE)
			output.write("\n");
		else
			output.write(value+"\n");
	}

	private static void writeDoubleRow(BufferedWriter output, String[] rowLabel, 
	                                   Matrix mat, int row, String delimiter) throws IOException {
		int nCols = mat.getNCols();
		for (int hdr=0; hdr < rowLabel.length; hdr++) {
		  output.write(quote(rowLabel[hdr], delimiter)+delimiter);
    }
		for (int col = 0; col < nCols-1; col++) {
			Double value = (Double)mat.getValue(row, col);
			if (value == null || Double.isNaN(value))
				output.write(delimiter);
			else
				output.write(value+delimiter);
		}
		Double value = (Double)mat.getValue(row, nCols-1);
		if (value == null || Double.isNaN(value))
			output.write("\n");
		else
			output.write(value+"\n");
	}
}
