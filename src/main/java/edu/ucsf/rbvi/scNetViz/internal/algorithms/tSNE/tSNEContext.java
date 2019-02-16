package edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE;

import java.util.List;

import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.BoundedDouble;

public class tSNEContext implements TSneConfiguration {
	public boolean cancelled = false;
	List<String> rowLabels = null;
	List<String> columnLabels = null;

	@Tunable(description="Initial Dimensions", 
	         longDescription="The number of dimensions to reduce the data set to before running "+
	                         "tSNE.  If the dimensionality of the data exceeds this number, "+
	                         "Principal Component Analysis (pca) will be used to do an initial "+
	                         "dimensionality reduction.  Settings this value to -1 ensures that "+
	                         "pca is not called.",
	         exampleStringValue="30",
	         groups={"t-SNE Parameters"}, gravity=66, format="#0")
	public int dimensions=10;

	@Tunable(description="Perplexity", 
	         longDescription="Perplexity is the balance between the local and global aspects of the data.",
	         exampleStringValue="20",
	         groups={"t-SNE Parameters"}, gravity=67)
	public double perplexity=20;

	@Tunable(description="Number of Iterations", 
	         longDescription="The number of iterations of the algorithm to perform",
	         exampleStringValue="5000",
	         groups={"t-SNE Parameters"}, gravity=68)
	public int iterations=5000;

	@Tunable(description="Use Barnes-Hut approximation", 
	         longDescription="The Barnes-Hut approximation is a way to reduce the computational complexity "+
	                         "of an algorithm by replacing a group of distant nodes with a single node at the "+
	                         "center of mass of all of those nodes",
	         exampleStringValue="false",
	         groups={"t-SNE Parameters"}, gravity=69)
	public boolean useBarnesHut=true;

	/*
	@Tunable(description="Theta value for Barnes-Hut", 
	         longDescription="The threshold value to activate Barnes-Hut.  This value reflects the accuracy "+
	                         "of the simulation.  If theta=0 then the approximation is never used",
	         exampleStringValue="0.5",
	         dependsOn="useBarnesHut=true", groups={"t-SNE Parameters"}, gravity=70)
	public BoundedDouble theta=new BoundedDouble(0.0, 0.9, 2.0, false, false);
	*/

	@Tunable(description="Theta value for Barnes-Hut", 
	         longDescription="The threshold value to activate Barnes-Hut.  This value reflects the accuracy "+
	                         "of the simulation.  If theta=0 then the approximation is never used",
	         exampleStringValue="0.001",
	         groups={"t-SNE Parameters"}, gravity=70)
	public BoundedDouble theta=new BoundedDouble(0.0, 0.001, 2.0, false, false);

	@Tunable(description="Log normalize the data", 
	         longDescription="Normalize the data by taking the log of each data point",
	         exampleStringValue="true",
	         groups={"t-SNE Parameters"}, gravity=76)
	public boolean logNormalize=true;

	@Tunable(description="Center and scale the data", 
	         longDescription="Center and scale the data before calculating the tSNE",
	         exampleStringValue="true",
	         groups={"t-SNE Parameters"}, gravity=76)
	public boolean centerAndScale=false;

	public tSNEContext(){}

	double[][] Xin = null;
	public double[][] getXin() {
		return Xin;
	}

	public void setXin(double[][] xin) {
		Xin = xin;
	}

	int outputDims = 3;
	public int getOutputDims() {
		return outputDims;
	}

	public void setOutputDims(int n) {
		outputDims = n;
	}

	public int getInitialDims() {
		return dimensions;
	}

	public void setInitialDims(int initial_dims) {
		dimensions = initial_dims;
	}

	public double getPerplexity() {
		return perplexity;
	}

	public void setPerplexity(double perplexity) {
		this.perplexity = perplexity;
	}

	public int getMaxIter() {
		return iterations;
	}

	public void setMaxIter(int max_iter) {
		iterations = max_iter;
	}

	public double getTheta() {
		return theta.getValue();
	}

	public void setTheta(double theta) {
		this.theta.setValue(theta);
	}

	boolean silent = false;
	public boolean silent() {
		return silent;
	}

	public void setSilent(boolean silent) {
		this.silent = silent;
	}

	public boolean printError() {
		return true;
	}

	public void setPrintError(boolean print_error) {
	}

	public int getXStartDim() {
		if (Xin != null && Xin[0] != null) return Xin[0].length;
		return 0;
	}

	public int getNrRows() {
		if (Xin != null) return Xin.length;
		return 0;
	}

	public boolean cancelled() {
		return cancelled;
	}

	@Override
	public boolean usePca() {
		return true;
	}

	@Override
	public void setUsePca(boolean use_pca) {}

	@Override
	public boolean centerAndScale() { return centerAndScale; }

	@Override
	public boolean logNormalize() { return logNormalize; }

	@Override
	public List<String> getRowLabels() { return rowLabels; }

	@Override
	public void setRowLabels(List<String> labels) { this.rowLabels = labels; }

	@Override
	public List<String> getColumnLabels() { return columnLabels; }

	@Override
	public void setColumnLabels(List<String> labels) { this.columnLabels = labels; }
}
