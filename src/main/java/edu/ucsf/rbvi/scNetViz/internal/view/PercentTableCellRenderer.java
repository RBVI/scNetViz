package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Component;
import java.text.NumberFormat;

import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;

public class PercentTableCellRenderer extends DefaultTableCellRenderer {
	NumberFormat formatter;

	public PercentTableCellRenderer() { super(); }

	public void setValue(Object value) {
		if (formatter == null) {
			formatter = NumberFormat.getPercentInstance();
		}
		setText((value == null || Double.isNaN((Double)value)) ? "" : formatter.format(value));
		setHorizontalAlignment(SwingConstants.RIGHT);
	}
}
