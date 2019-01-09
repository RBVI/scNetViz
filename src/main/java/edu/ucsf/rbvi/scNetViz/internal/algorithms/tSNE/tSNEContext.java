package edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE;

import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.BoundedDouble;


public class tSNEContext implements TSneConfiguration {
	public boolean cancelled = false;

	@Tunable(description="Initial Dimensions", 
	         longDescription="The number of dimensions to reduce the data set to before running "+
	                         "tSNE.  If the dimensionality of the data exceeds this number, "+
	                         "Principal Component Analysis (pca) will be used to do an initial "+
	                         "dimensionality reduction.  Settings this value to -1 ensures that "+
	                         "pca is not called.",
	         exampleStringValue="30",
	         groups={"t-SNE Parameters"}, gravity=66, format="#0")
	public int dimensions=-1;

	@Tunable(description="Perplexity", 
	         longDescription="Perplexity is the balance between the local and global aspects of the data.",
	         exampleStringValue="20",
	         groups={"t-SNE Parameters"}, gravity=67)
	public double perplexity=20;

	@Tunable(description="Number of Iterations", 
	         longDescription="The number of iterations of the algorithm to perform",
	         exampleStringValue="2000",
	         groups={"t-SNE Parameters"}, gravity=68)
	public int iterations=2000;

	@Tunable(description="Use Barnes-Hut approximation", 
	         longDescription="The Barnes-Hut approximation is a way to reduce the computational complexity "+
	                         "of an algorithm by replacing a group of distant nodes with a single node at the "+
	                         "center of mass of all of those nodes",
	         exampleStringValue="false",
	         groups={"t-SNE Parameters"}, gravity=69)
	public boolean useBarnesHut=false;

	@Tunable(description="Theta value for Barnes-Hut", 
	         longDescription="The threshold value to activate Barnes-Hut.  This value reflects the accuracy "+
	                         "of the simulation.  If theta=0 then the approximation is never used",
	         exampleStringValue="0.9",
	         dependsOn="useBarnesHut=true", groups={"t-SNE Parameters"}, gravity=70)
	public BoundedDouble theta=new BoundedDouble(0.0, 0.9, 2.0, false, false);

	public tSNEContext(){}

	double[][] Xin = null;
	public double[][] getXin() {
		return Xin;
	}

	public void setXin(double[][] xin) {
		Xin = xin;
	}

	int outputDims = 2;
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
		return false;
	}

	@Override
	public void setUsePca(boolean use_pca) {}
}
