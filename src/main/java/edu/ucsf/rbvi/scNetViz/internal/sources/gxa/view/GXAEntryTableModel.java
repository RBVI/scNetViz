package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.view;

import java.util.List;

import javax.swing.table.AbstractTableModel;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;

public class GXAEntryTableModel extends AbstractTableModel {
	List<Metadata> entries = null;

	static String[] columnNames = {"Accession", "Loaded", "Experiment", "Assays", "Comparisons", "Organisms", "Experimental Variables"};

	public GXAEntryTableModel (List<Metadata> entries) {
		super();
		this.entries = entries;
	}

	@Override
	public int getColumnCount() { return columnNames.length; }

	@Override
	public String getColumnName(int column) {
		return columnNames[column];
	}

	@Override
	public int getRowCount() { 
		return entries.size(); 
	}

	@Override
	public Class getColumnClass(int column) {
		switch (column) {
			case 0:
				return String.class;
			case 1:
				return String.class;
			case 2:
				return String.class;
			case 3:
				return Integer.class;
			case 4:
				return Integer.class;
			case 5:
				return String.class;
			case 6:
				return String.class;
		}
		return String.class;
	}

	@Override
	public Object getValueAt(int row, int column) {
		GXAMetadata entry = (GXAMetadata)entries.get(row);
		switch (column) {
			case 0:
				return entry.get(GXAMetadata.ACCESSION);
			case 1:
				return entry.get(GXAMetadata.DATE);
			case 2:
				return entry.get(GXAMetadata.DESCRIPTION);
			case 3:
				return new Integer(((Long)entry.get(GXAMetadata.ASSAYS)).intValue());
			case 4:
				return new Integer(((Long)entry.get(GXAMetadata.CONTRASTS)).intValue());
			case 5:
				return entry.get(GXAMetadata.SPECIES);
			case 6:
				return nlList(entry.get(GXAMetadata.FACTORS));
		}
		return null;
	}

	private String nlList(Object list) {
		String ret = null;
		if (!(list instanceof List)) {
			return null;
		}
		List<String> strList = (List<String>)list;
		for (String l: strList) {
			if (ret == null)
				ret = "<html>"+l;
			else
				ret += "<br/>"+l;
		}
		return ret+"</html>";
	}
}
