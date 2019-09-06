package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskObserver;

public class Species implements Comparable<Species> {
	private static List<Species> allSpecies;
	private int taxon_id;
	private String compactName;
	private String officialName;

	public Species(int tax, String name, String oName) {
		this.taxon_id = tax;
		this.compactName = name;
		this.officialName = oName;
	}

	public String toString() { return compactName; }
	public String getName() { return compactName; }
	public String getOfficialName() { return officialName; }
	public int getTaxId() { return taxon_id; }

	@Override
	public int compareTo(Species t) {
		if (t.toString() == null) return 1;
		return this.toString().compareTo(t.toString());
	}

	public static List<Species> getSpecies() { return allSpecies; }

	public static void loadSpecies(final ScNVManager manager) {
		LoadSpecies ls = new LoadSpecies(manager);
		Thread t = new Thread(ls);
		t.start();
	}

	static class LoadSpecies implements Runnable, TaskObserver {
		final ScNVManager manager;
		LoadSpecies(final ScNVManager manager) { 
			this.manager = manager; 
		}

		@Override
		public void run() {
			while (!manager.haveCommand("string", "list species")) {
				try {
					Thread.sleep(500);
				} catch (Exception e) {}
			}
			if (manager.haveCommand("string", "list species"))
				manager.executeCommand("string", "list species", null, this, true);
		}

		@Override
		public void allFinished(FinishStatus status) {}

		@Override
		public void taskFinished(ObservableTask task) {
			allSpecies = new ArrayList<Species>();
			List<Map<String,String>> list = task.getResults(List.class);
			if (list == null) {
				System.out.println("Ack!  No species");
				return;
			}
			for (Map<String,String> specie: list) {
				int taxId = Integer.parseInt(specie.get("taxonomyId"));
				Species sp = new Species(taxId, specie.get("abbreviatedName"), specie.get("scientificName"));
				allSpecies.add(sp);
			}
		}
	}
}
