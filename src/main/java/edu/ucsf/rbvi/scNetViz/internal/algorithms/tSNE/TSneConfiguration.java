package edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE;

import java.util.List;

public interface TSneConfiguration {

	double[][] getXin();

	void setXin(double[][] xin);

	int getOutputDims();

	void setOutputDims(int n);

	int getInitialDims();

	void setInitialDims(int initial_dims);

	double getPerplexity();

	void setPerplexity(double perplexity);

	int getMaxIter();

	void setMaxIter(int max_iter);

	boolean usePca();

	void setUsePca(boolean use_pca);

	double getTheta();

	void setTheta(double theta);

	boolean silent();

	void setSilent(boolean silent);

	boolean printError();

	void setPrintError(boolean print_error);

	int getXStartDim();

	int getNrRows();

	boolean cancelled();

	boolean centerAndScale();

	boolean logNormalize();

	void setRowLabels(List<String> labels);

	List<String> getRowLabels();

	void setColumnLabels(List<String> labels);

	List<String> getColumnLabels();
}
