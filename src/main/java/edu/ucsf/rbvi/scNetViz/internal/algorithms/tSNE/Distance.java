package edu.ucsf.rbvi.scNetViz.internal.algorithms.tSNE;

public interface Distance {
	double distance(DataPoint d1, DataPoint d2);
}
