package edu.ucsf.rbvi.scNetViz.internal.tasks;

import org.cytoscape.work.Tunable;

public class AdvancedRemoteParameters {

	@Tunable (description="Minimum number of genes/cell",
	          groups={"Advanced preprocessing parameters"},params="displayState=collapsed",
	          tooltip="Filter cells that don't have at least this number of genes")
	public int min_genes=100;

	@Tunable (description="Minimum number of cells expressed by a gene",
	          groups={"Advanced preprocessing parameters"},params="displayState=collapsed",
	          tooltip="Filter genes not expressed in at least this number of cells")
	public int min_cells=1;

	@Tunable (description="Normalize",
	          groups={"Advanced preprocessing parameters"},params="displayState=collapsed",
	          tooltip="Normalize the data based on the total TPMs")
	public boolean normalize=true;

	@Tunable (description="Log transform",
	          groups={"Advanced preprocessing parameters"},params="displayState=collapsed",
	          tooltip="Log transform the data")
	public boolean log1p=true;

	@Tunable (description="Highly variable genes",
	          groups={"Advanced preprocessing parameters"},params="displayState=collapsed",
	          tooltip="Restrict the data to only highly variable genes")
	public boolean hvg=true;

	@Tunable (description="Scale the final matrix",
	          groups={"Advanced preprocessing parameters"},params="displayState=collapsed",
	          tooltip="Scale the data final matrix after all preprocessing steps")
	public boolean scale=true;

	public String getArgs() {
		StringBuilder builder = new StringBuilder();
		builder.append("min_genes="+min_genes);
		builder.append("&min_cells="+min_cells);
		builder.append("&normalize="+normalize);
		builder.append("&log1p="+log1p);
		builder.append("&hvg="+hvg);
		builder.append("&scale="+scale);
		return builder.toString();
	}
}
