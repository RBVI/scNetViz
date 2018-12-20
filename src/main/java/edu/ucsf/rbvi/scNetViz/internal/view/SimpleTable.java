package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import java.util.List;

import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

import edu.ucsf.rbvi.scNetViz.internal.api.PercentDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.MyDouble;
import edu.ucsf.rbvi.scNetViz.internal.api.PValueDouble;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class SimpleTable extends JTable {
	final ScNVManager manager;
	final SortableTableModel tableModel;
	final SimpleTable thisComponent;

	static Color alternateColor = new Color(234,255,234);

	public SimpleTable (final ScNVManager manager, final SortableTableModel tableModel) {
		super(tableModel);
		this.manager = manager;
		this.tableModel = tableModel;
		this.thisComponent = this;

		this.setAutoCreateRowSorter(true);
		this.setAutoCreateColumnsFromModel(true);
		// this.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
		this.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		int nColumns = tableModel.getColumnCount();
		TableColumnModel columnModel = getColumnModel();
		columnModel.getColumn(0).setPreferredWidth(100);
		for (int col = 1; col < nColumns; col++) {
			columnModel.getColumn(col).setPreferredWidth(100);
		}

		setDefaultRenderer(MyDouble.class, new DoubleTableCellRenderer());
		setDefaultRenderer(PercentDouble.class, new PercentTableCellRenderer());
		setDefaultRenderer(PValueDouble.class, new PValueTableCellRenderer());

		JTableHeader header = getTableHeader();
		header.setDefaultRenderer(new HeaderRenderer(this));
		header.setFont(new Font("SansSerif", Font.BOLD, 10));

		setRowSorter(new NaNTableRowSorter(tableModel));

		addMouseListener(new MouseAdapter() {
			public void mousePressed(MouseEvent mouseEvent) {
				JTable table = (JTable) mouseEvent.getSource();
				Point point = mouseEvent.getPoint();
				int row = table.rowAtPoint(point);
				int col = table.columnAtPoint(point);
				if (col != 0) return;
				if (mouseEvent.getClickCount() == 2 && table.getSelectedRow() != -1) {
					tableModel.sortColumns(row);
				}
			}
		});

		this.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		ListSelectionModel selectionModel = getSelectionModel();

		selectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				int row = thisComponent.getSelectedRow();
				if (row >= 0)
					tableModel.setSelectedRow(thisComponent.convertRowIndexToModel(row));
			}
		});

		doLayout();
	}

	@Override
	public void tableChanged(TableModelEvent e) {
		super.tableChanged(e);
		if (tableModel == null) return;
		TableColumnModel columnModel = getColumnModel();
		columnModel.getColumn(0).setPreferredWidth(100);
		for (int col = 1; col < tableModel.getColumnCount(); col++) {
			columnModel.getColumn(col).setPreferredWidth(100);
		}
	}

	@Override
	public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
		Component returnComp = super.prepareRenderer(renderer, row, column);

		if (!returnComp.getBackground().equals(getSelectionBackground())) {
			Color bg = (row % 2 == 0 ? alternateColor : Color.WHITE);
			returnComp.setBackground(bg);
			returnComp.setFont(new Font("SansSerif", Font.PLAIN, 10));
			bg = null;
		}
		return returnComp;
	}
}
