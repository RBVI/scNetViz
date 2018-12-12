package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Component;
import java.text.DecimalFormat;

import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;

public class DoubleTableCellRenderer extends DefaultTableCellRenderer {
	DecimalFormat formatter;

	public DoubleTableCellRenderer() { super(); }

	public void setValue(Object value) {
		if (formatter == null) {
			formatter = new DecimalFormat("###,###.000");
		}
		setText((value == null || Double.isNaN((Double)value)) ? "" : formatter.format(value));
		setHorizontalAlignment(SwingConstants.RIGHT);
	}
}
