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
				double[][] mat = ((DoubleMatrix)matrix).getDoubleMatrix(Double.NaN);
				for (int row = 0; row < rowLabels.size(); row++) {
					writeRow(output, rowLabels.get(row), mat[row], delimiter);
				}
			} else if (matrix instanceof IntegerMatrix) {
				int[][] mat = ((IntegerMatrix)matrix).getIntegerMatrix(Integer.MIN_VALUE);
				for (int row = 0; row < rowLabels.size(); row++) {
					writeRow(output, rowLabels.get(row), mat[row], delimiter);
				}
			} else if (matrix instanceof StringMatrix) {
				String[][] mat = ((StringMatrix)matrix).getStringMatrix();
				for (int row = 0; row < rowLabels.size(); row++) {
					writeRow(output, rowLabels.get(row), mat[row], delimiter);
				}
			}
		
			// Write it all out
		} catch (Exception fnf) {
		}
	}

	private static String quote(String str, String delimiter) {
		if (delimiter.equals("\t"))
			return str;
		return "\""+str+"\"";
	}

	private static void writeRow(BufferedWriter output, String rowLabel, String[] row, String delimiter) throws IOException {
		output.write(quote(rowLabel, delimiter)+delimiter);
		for (int i = 0; i < row.length-1; i++) {
			output.write(quote(row[i], delimiter)+delimiter);
		}
		output.write(quote(row[row.length-1], delimiter)+"\n");
	}

	private static void writeRow(BufferedWriter output, String rowLabel, int[] row, String delimiter) throws IOException {
		output.write(quote(rowLabel, delimiter)+delimiter);
		for (int i = 0; i < row.length-1; i++) {
			if (row[i] == Integer.MIN_VALUE)
				output.write(delimiter);
			else
				output.write(row[i]+delimiter);
		}
		output.write(row[row.length-1]+"\n");
	}

	private static void writeRow(BufferedWriter output, String rowLabel, double[] row, String delimiter) throws IOException {
		output.write(quote(rowLabel, delimiter)+delimiter);
		for (int i = 0; i < row.length-1; i++) {
			if (Double.isNaN(row[i]))
				output.write(delimiter);
			else
				output.write(row[i]+delimiter);
		}
		output.write(row[row.length-1]+"\n");
	}
}
