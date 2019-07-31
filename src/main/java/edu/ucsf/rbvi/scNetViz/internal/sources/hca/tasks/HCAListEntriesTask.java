package edu.ucsf.rbvi.scNetViz.internal.sources.hca.tasks;

import java.util.Arrays;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.json.JSONResult;

import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.hca.HCASource;

public class HCAListEntriesTask extends AbstractTask implements ObservableTask {
	final ScNVManager scManager;
	final HCASource hcaSource;

	@Tunable(description="Species to include in list", context="nogui", required=false,
	         longDescription="Species name.  This should be the actual "+
                           "taxonomic name (e.g. homo sapiens, not human)",
           exampleStringValue="homo sapiens")
	public String species = null;

	@Tunable(description="Minimum number of assays to include", context="nogui", required=false)
	public int assays = 0;

	public HCAListEntriesTask(final ScNVManager scManager, final HCASource hcaSource) {
		super();
		this.scManager = scManager;
		this.hcaSource = hcaSource;
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		taskMonitor.setTitle("List Human Cell Atlas Entries");
		taskMonitor.setStatusMessage("List Human Cell Atlas");
	}

	@Override
	public List<Class<?>> getResultClasses() {
		return Arrays.asList(String.class, JSONResult.class, List.class);
	}

	@Override
	public <R> R getResults(Class<? extends R> type) {
		List<Metadata> metadata = hcaSource.getMetadata();
		if (type.equals(JSONResult.class)) {
			return (R)(new JSONEnumeration(metadata));
		} else if (type.equals(List.class)) {
			return (R)metadata;
		} else {
			// Default to String.class
			StringBuilder builder = new StringBuilder();
			builder.append("Available experiments: \n");
			for (Metadata meta: metadata) {
				if (species != null &&
						!species.equalsIgnoreCase(meta.get(Metadata.SPECIES).toString()))
					continue;
				if (Integer.parseInt(meta.get(HCAMetadata.ASSAYS).toString()) < assays)
					continue;

				builder.append(meta.get(Metadata.ACCESSION)+": ");
				builder.append(meta.get(Metadata.DESCRIPTION));
				builder.append(" from "+meta.get(Metadata.SPECIES)+" ");
				builder.append(" with "+meta.get(HCAMetadata.ASSAYS)+" assays");
				builder.append("\n");
			}
			return (R)builder.toString();
		}
	}

	class JSONEnumeration implements JSONResult {
		final List<Metadata> metadata;
		JSONEnumeration(final List<Metadata> metadata) {
			this.metadata = metadata;
		}
	
		public String getJSON() {
			StringBuilder builder = new StringBuilder();
			builder.append("[");
			for (Metadata meta: metadata) {
				if (species != null &&
						!species.equalsIgnoreCase(meta.get(Metadata.SPECIES).toString()))
					continue;
				if (Integer.parseInt(meta.get(HCAMetadata.ASSAYS).toString()) < assays)
					continue;
				builder.append(((HCAMetadata)meta).toJSON()+",");
			}
			return builder.substring(0, builder.length()-1)+"]";
		}
	}

	@ProvidesTitle
	public String getTitle() {return "List Human Cell Atlas Entries";}
}
