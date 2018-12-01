package edu.ucsf.rbvi.scNetViz.internal.utils;

import java.util.Arrays;
import java.util.Comparator;

public class MatrixUtils {

	public static Integer[] indexSort(double[] tData, int nVals) {
    Integer[] index = new Integer[nVals];
    for (int i = 0; i < nVals; i++) index[i] = i;
    IndexComparator iCompare = new IndexComparator(tData);
    Arrays.sort(index, iCompare);
    return index;
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

    public IndexComparator(String[] data) { this.stringData = data; }

    public IndexComparator(double[] data) { this.data = data; }

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

        if (data[o1] < data[o2]) return -1;
        if (data[o1] > data[o2]) return 1;
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
