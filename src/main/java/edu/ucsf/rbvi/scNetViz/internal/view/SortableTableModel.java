package edu.ucsf.rbvi.scNetViz.internal.view;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import javax.swing.ListSelectionModel;
import javax.swing.event.TableModelEvent;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

public abstract class SortableTableModel extends AbstractTableModel {
	protected Integer[] columnIndex = null;
	protected int hdrCols;
	private int lastRow = -1;
	private int lastDirection = 0;

	public SortableTableModel(int hdrCols) {
		super();
		this.hdrCols = hdrCols;
	}

	abstract public int getSelectedRow();
	abstract public void setSelectedRow(int selectedRow);

	public void sortColumns(int row) {
		System.out.println("sorting row: "+row);
		System.out.println("hdrCols = "+hdrCols);
		System.out.println("getColumnCount()-hdrCols = "+(getColumnCount()-hdrCols));
		Integer tmpIndex[] = new Integer[getColumnCount()-hdrCols];
		columnIndex = null;
		for (int i = 0; i < (getColumnCount()-hdrCols); i++) 
			tmpIndex[i] = i;
		IndexSorter sorter = new IndexSorter(this, row, hdrCols);
		Arrays.sort(tmpIndex, sorter);
		columnIndex = tmpIndex;
		// for (int i = 0; i < getColumnCount()-hdrCols; i++) {
		// 	System.out.println("i = "+i+" columnIndex[i] = "+columnIndex[i]);
		// 	System.out.println(getColumnName(i));
		// 	System.out.println(getValueAt(row,i));
		// }

		fireTableChanged(new TableModelEvent(this, TableModelEvent.HEADER_ROW));
	}

	public class IndexSorter implements Comparator<Integer> {
		final TableModel tableModel;
		int row;
		int hdrCols;
		Class<?> rowType;

		public IndexSorter(final TableModel tableModel, int row, int hdrCols) {
			this.tableModel = tableModel;
			this.row = row;
			if (row != lastRow) {
				lastRow = row;
				lastDirection = -1; // High to low
			} else {
				lastDirection = lastDirection * (-1);
			}

			this.hdrCols = hdrCols;
			rowType = tableModel.getColumnClass(hdrCols);
		}

		public int compare(Integer o1, Integer o2) {
			// if (o1 < hdrCols || o2 < hdrCols) 
			// 	return o1.compareTo(o2);
			// if (true) return 0;
			// Note: these are getting the data from the JTable, *not* the Matrix
			Object v1 = tableModel.getValueAt(row, o1+hdrCols);
			Object v2 = tableModel.getValueAt(row, o2+hdrCols);
			if (lastDirection < 0) {
				Object tmp = v2;
				v2 = v1; v1 = tmp;
			}
			if (v1 == null && v2 == null) return 0;
			if (v1 == null) return -1;
			if (v2 == null) return 1;
			if (rowType.equals(String.class)) {
				return ((String)v1).compareTo((String)v2);
			} else if (rowType.equals(Integer.class)) {
				return ((Integer)v1).compareTo((Integer)v2);
			} else if (rowType.equals(Double.class)) {
				return ((Double)v1).compareTo((Double)v2);
			}
			return v1.toString().compareTo(v2.toString());
		}

	}
}
