package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;

public class MatrixUtils {

	public static String CONTROL_PREFIX = "ERCC";

	public static BitSet findControls(final List<String[]> rowLabels, final String prefix, int hdr) {
		BitSet bitSet = new BitSet(rowLabels.size());
		for (int index = 0; index < rowLabels.size(); index++) {
			if (rowLabels.get(index)[hdr].startsWith(prefix))
				bitSet.set(index);
		}
		return bitSet;
	}

	public static List<String> filterLabels(final List<String[]> rowLabels, BitSet excludeRows, int hdr) {
		List<String> newList = new ArrayList<String>(rowLabels.size() - excludeRows.cardinality());
		for (int row = 0; row < rowLabels.size(); row++) {
			if (excludeRows.get(row))
				continue;
			newList.add(rowLabels.get(row)[hdr]);
		}
		return newList;
	}

	public static Integer[] indexSort(double[] tData, int nVals, boolean abs) {
    Integer[] index = new Integer[nVals];
    for (int i = 0; i < nVals; i++) index[i] = i;
    IndexComparator iCompare = new IndexComparator(tData, abs);
    Arrays.sort(index, iCompare);
    return index;
	}

	public static Integer[] indexSort(double[] tData, int nVals) {
		return indexSort(tData, nVals, false);
  }

	public static Integer[] indexSort(int[] tData, int nVals) {
    Integer[] index = new Integer[nVals];
    for (int i = 0; i < nVals; i++) index[i] = i;
    IndexComparator iCompare = new IndexComparator(tData);
    Arrays.sort(index, iCompare);
    return index;
  }

	public static Integer[] indexSort(String[] tData, int nVals) {
    Integer[] index = new Integer[nVals];
    for (int i = 0; i < nVals; i++) index[i] = i;
    IndexComparator iCompare = new IndexComparator(tData);
    Arrays.sort(index, iCompare);
    return index;
  }


	private static class IndexComparator implements Comparator<Integer> {
    double[] data = null;
    int[] intData = null;
    String[] stringData = null;
		boolean abs = false;

    public IndexComparator(String[] data) { this.stringData = data; }

    public IndexComparator(double[] data) { this.data = data; }

    public IndexComparator(double[] data, boolean abs) { 
			this.data = data; 
			this.abs = abs;
		}

    public IndexComparator(int[] data) { this.intData = data; }

    public int compare(Integer o1, Integer o2) {
      if (data != null) {
				// NaN handling
				if (Double.isNaN(data[o1]) && Double.isNaN(data[o2]))
					return 0;
				if (Double.isNaN(data[o1]))
					return -1;
				if (Double.isNaN(data[o2]))
					return 1;

				if (abs) {
        	if (Math.abs(data[o1]) < Math.abs(data[o2])) return -1;
        	if (Math.abs(data[o1]) > Math.abs(data[o2])) return 1;
				} else {
        	if (data[o1] < data[o2]) return -1;
        	if (data[o1] > data[o2]) return 1;
				}
        return 0;
      } else if (intData != null) {
        if (intData[o1] < intData[o2]) return -1;
        if (intData[o1] > intData[o2]) return 1;
        return 0;
      } else if (stringData != null) {
        return stringData[o1].compareTo(stringData[o2]);
      }
      return 0;
    }
  }
}
