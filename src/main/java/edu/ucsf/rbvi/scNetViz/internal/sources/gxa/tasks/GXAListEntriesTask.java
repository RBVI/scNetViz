package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.tasks;

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
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;

public class GXAListEntriesTask extends AbstractTask implements ObservableTask {
	final ScNVManager scManager;
	final GXASource gxaSource;

	@Tunable(description="Species to include in list")
	public String species;

	@Tunable(description="Minimum number of assays to include")
	public int assays = 0;

	public GXAListEntriesTask(final ScNVManager scManager, final GXASource gxaSource) {
		super();
		this.scManager = scManager;
		this.gxaSource = gxaSource;
	}

	@Override
	public void run(TaskMonitor taskMonitor) {
		taskMonitor.setTitle("List Single Cell Expression Atlas Entries");
		taskMonitor.setStatusMessage("List Single Cell Expression Atlas");
	}

	@Override
	public List<Class<?>> getResultClasses() {
		return Arrays.asList(String.class, JSONResult.class, List.class);
	}

	@Override
	public <R> R getResults(Class<? extends R> type) {
		List<Metadata> metadata = gxaSource.getMetadata();
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
				if (Integer.parseInt(meta.get(GXAMetadata.ASSAYS).toString()) < assays)
					continue;

				builder.append(meta.get(Metadata.ACCESSION)+": ");
				builder.append(meta.get(Metadata.DESCRIPTION));
				builder.append(" from "+meta.get(Metadata.SPECIES)+" ");
				builder.append(" with "+meta.get(GXAMetadata.ASSAYS)+" assays");
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
				builder.append(((GXAMetadata)meta).toJSON()+",");
			}
			builder.append("]");
			return builder.toString();
		}
	}

	@ProvidesTitle
	public String getTitle() {return "List Single Cell Expression Atlas Entries";}
}
