package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Component;
import java.text.DecimalFormat;

import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;

public class PValueTableCellRenderer extends DefaultTableCellRenderer {
	DecimalFormat formatter;

	public PValueTableCellRenderer() { super(); }

	public void setValue(Object value) {
		if (formatter == null) {
			formatter = new DecimalFormat("0.00E0");
		}
		setText((value == null || Double.isNaN((Double)value)) ? "" : formatter.format(value));
		setHorizontalAlignment(SwingConstants.RIGHT);
	}
}
