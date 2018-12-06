package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.StringMatrix;

public class CSVWriter {
	public static void writeCSV(File file, Matrix matrix, String delimiter) {
		try {
			BufferedWriter output = new BufferedWriter(new FileWriter(file));
			List<String> rowLabels = matrix.getRowLabels();
			// Get the column labels
			List<String> colLabels = matrix.getColLabels();

			for (int i = 0; i < colLabels.size()-1; i++) {
				output.write(quote(colLabels.get(i), delimiter)+delimiter);
			}
			output.write(quote(colLabels.get(colLabels.size()-1), delimiter)+"\n");

			// Get the matrix
			if (matrix instanceof DoubleMatrix) {
				for (int row = 0; row < rowLabels.size(); row++) {
					writeDoubleRow(output, rowLabels.get(row), (DoubleMatrix)matrix, row, delimiter);
				}
			} else if (matrix instanceof IntegerMatrix) {
				for (int row = 0; row < rowLabels.size(); row++) {
					writeIntegerRow(output, rowLabels.get(row), (IntegerMatrix)matrix, row, delimiter);
				}
			} else if (matrix instanceof StringMatrix) {
				for (int row = 0; row < rowLabels.size(); row++) {
					writeStringRow(output, rowLabels.get(row), (StringMatrix)matrix, row, delimiter);
				}
			}
		
			// Write it all out
		} catch (Exception fnf) {
			fnf.printStackTrace();
		}
	}

	private static String quote(String str, String delimiter) {
		if (delimiter.equals("\t"))
			return str;
		return "\""+str+"\"";
	}

	private static void writeStringRow(BufferedWriter output, String rowLabel, 
	                                   StringMatrix mat, int row, String delimiter) throws IOException {
		int nCols = mat.getNCols();
		output.write(quote(rowLabel, delimiter)+delimiter);
		for (int col = 0; col < nCols-1; col++) {
			output.write(quote(mat.getValue(row, col), delimiter)+delimiter);
		}
		output.write(quote(mat.getValue(row,nCols-1), delimiter)+"\n");
	}

	private static void writeIntegerRow(BufferedWriter output, String rowLabel, 
	                                    IntegerMatrix mat, int row, String delimiter) throws IOException {
		int nCols = mat.getNCols();
		output.write(quote(rowLabel, delimiter)+delimiter);
		for (int col = 0; col < nCols-1; col++) {
			int value = mat.getIntegerValue(row, col);
			if (value == Integer.MIN_VALUE)
				output.write(delimiter);
			else
				output.write(value+delimiter);
		}
		output.write(mat.getIntegerValue(row,nCols-1)+"\n");
	}

	private static void writeDoubleRow(BufferedWriter output, String rowLabel, 
	                                   DoubleMatrix mat, int row, String delimiter) throws IOException {
		int nCols = mat.getNCols();
		output.write(quote(rowLabel, delimiter)+delimiter);
		for (int col = 0; col < nCols-1; col++) {
			double value = mat.getDoubleValue(row, col);
			if (Double.isNaN(value))
				output.write(delimiter);
			else
				output.write(value+delimiter);
		}
		output.write(mat.getDoubleValue(row,nCols-1)+"\n");
	}
}
