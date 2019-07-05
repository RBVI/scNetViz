package edu.ucsf.rbvi.scNetViz.internal.sources.hca.view;

import java.util.ArrayList;
import java.util.List;

import javax.swing.table.AbstractTableModel;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCASource;

public class HCAEntryTableModel extends AbstractTableModel {
	List<Metadata> entries = null;

	static String[] columnNames = {"Accession", "Type", "Description", "Cells", "Organisms", "Tissues"};

	public HCAEntryTableModel (List<Metadata> entries) {
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
				return String.class;
			case 5:
				return String.class;
		}
		return String.class;
	}

	@Override
	public Object getValueAt(int row, int column) {
		HCAMetadata entry = (HCAMetadata)entries.get(row);
		switch (column) {
			case 0:
				return entry.get(HCAMetadata.ACCESSION);
			case 1:
				return entry.get(HCAMetadata.TYPE);
			case 2:
				return entry.get(HCAMetadata.DESCRIPTION);
			case 3:
				return new Integer(((Long)entry.get(HCAMetadata.ASSAYS)).intValue());
			case 4:
				return entry.get(HCAMetadata.SPECIES);
			case 5:
				return entry.get(HCAMetadata.ORGANS);
		}
		return null;
	}

	public List<Integer> search(String searchText) {
		String lcString = searchText.toLowerCase();
		List<Integer> results = new ArrayList<>();
		for (int i = 0; i < entries.size(); i++) {
			Metadata entry = entries.get(i);
			if (contains(entry, HCAMetadata.ACCESSION, searchText) ||
			    contains(entry, HCAMetadata.DESCRIPTION, searchText) ||
					contains(entry, HCAMetadata.SPECIES, searchText))
				results.add(i);
		}
		return results;
	}

	public List<Integer> searchRegex(String searchText) {
		// Search all accessions, organisms, and descriptions
		List<Integer> results = new ArrayList<>();
		for (int i = 0; i < entries.size(); i++) {
			Metadata entry = entries.get(i);
			if (matches(entry, HCAMetadata.ACCESSION, searchText) ||
			    matches(entry, HCAMetadata.DESCRIPTION, searchText) ||
					matches(entry, HCAMetadata.SPECIES, searchText))
				results.add(i);
		}
		return results;
	}

	private boolean matches(Metadata entry, String key, String text) {
		return ((String)entry.get(key)).matches(text);
	}

	private boolean contains(Metadata entry, String key, String text) {
		String str = ((String)entry.get(key)).toLowerCase();
		return str.contains(text);
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
