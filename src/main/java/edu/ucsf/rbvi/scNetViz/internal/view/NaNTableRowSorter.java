package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Component;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Comparator;

import javax.swing.SwingConstants;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import edu.ucsf.rbvi.scNetViz.internal.api.MyDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.PercentDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.PValueDouble;

public class NaNTableRowSorter extends TableRowSorter<TableModel> {
	DecimalFormat formatter;
	TableModel myTableModel;
	Boolean[] sortOrder;
	int currentSortColumn = -1;
	Comparator<?> ascendingComparator;
	Comparator<?> descenndingComparator;

	public NaNTableRowSorter(TableModel tableModel) {
		super(tableModel);
		myTableModel = tableModel;
		sortOrder = new Boolean[tableModel.getColumnCount()];
		Arrays.fill(sortOrder, null);
	}

	@Override
	public void toggleSortOrder(int column) {
		if (sortOrder[column] == null)
			sortOrder[column] = true;
		else if (sortOrder[column]) {
			sortOrder[column] = false;
		} else {
			sortOrder[column] = true;
		}
		super.toggleSortOrder(column);
	}

	@Override
	public Comparator<?> getComparator(int column) {
		currentSortColumn = column;
		Class<?> columnClass = myTableModel.getColumnClass(column);
		if(columnClass.equals(MyDouble.class) || columnClass.equals(PercentDouble.class) ||
		   columnClass.equals(PValueDouble.class)) {
			return new NaNComparator(sortOrder[column]);
		}
		return super.getComparator(column);
	}

	public class NaNComparator implements Comparator {
		boolean ascending;

		public NaNComparator(boolean ascending) {
			this.ascending = ascending;
		}

		@Override
		public int compare(Object o1, Object o2) {
			if (o1 instanceof String && o2 instanceof String) {
				String s1 = (String)o1;
				String s2 = (String)o2;

				if (s1.compareTo(s2) == 0)
					return 0;

				if (s1.length() == 0) {
					if (ascending) return 1;
					return -1;
				}
				if (s2.length() == 0) {
					if (ascending) return -1;
					return 1;
				}

				Double d1 = Double.parseDouble(s1);
				Double d2 = Double.parseDouble(s2);

				return d1.compareTo(d2);
			}

			Double d1 = (Double)o1;
			Double d2 = (Double)o2;

			return Double.compare(d1,d2);
		}

	}
}
