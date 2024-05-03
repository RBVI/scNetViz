package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ContainsTunables;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.json.JSONResult;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.ExperimentFrame;
import edu.ucsf.rbvi.scNetViz.internal.view.ViewUtils;


// TODO: Consider random subsampling?
// TODO: Expose filter criteria
abstract public class AbstractEmbeddingTask extends AbstractTask implements ObservableTask{

	protected Map<String,double[]> embedding;
	final ScNVManager manager;
		
	public AbstractEmbeddingTask(final ScNVManager manager) {
		this.manager = manager;
	}

	public void cancel() {
		cancelled = true;
	}

	public Map<String,double[]> getResults() { return embedding; }

	protected Map<String,double[]> createEmbedding(List<String> lines) {
		embedding = new LinkedHashMap<String, double[]>();
		for (String line: lines) {
			String[] tokens = line.split(",");
			if (tokens.length <= 2) continue;
			double[] coords = new double[2];
			coords[0] = Double.valueOf(tokens[tokens.length-2]);
			coords[1] = Double.valueOf(tokens[tokens.length-1]);
			String cell = tokens[0];
			for (int i = 1; i < tokens.length-2; i++) {
				cell += ","+tokens[i];
			}
			embedding.put(cell, coords);
		}
		return embedding;
	}

	public void scale(Map<String, double[]> values) {
		double xMax = Double.MIN_VALUE;
		double xMin = Double.MAX_VALUE;
		double yMax = Double.MIN_VALUE;
		double yMin = Double.MAX_VALUE;
		for (String cell: values.keySet()) {
			double[] coords = values.get(cell);
			xMin = Math.min(xMin, coords[0]);
			xMax = Math.max(xMax, coords[0]);
			yMin = Math.min(yMin, coords[1]);
			yMax = Math.max(yMax, coords[1]);
		}
		double xScale = xMax - xMin;
		double yScale = yMax - yMin;

		for (String cell: values.keySet()) {
			double[] coords = values.get(cell);
			coords[0] = (coords[0]-xMin)/xScale;
			coords[1] = (coords[1]-yMin)/yScale;
		}
	}

	protected void createCache(MatrixMarket mmtx, Experiment exp) {
		String source = exp.getSource().getName();
		String accession = exp.getMetadata().get(Metadata.ACCESSION).toString();
		mmtx.createCache(source, accession);
		while (!mmtx.hasCache()) {
			try {
				Thread.sleep(500);
			} catch (Exception e) {
				break;
			}
		}
	}

  protected void updateCellPlot(Experiment exp, String plotType) {
    exp.setTSNE(embedding);
    exp.setPlotType(plotType);
	  ExperimentFrame frame = manager.getExperimentFrame(exp);
    frame.getTPMTab().updateCellPlot();
  }

	@Override
	public List<Class<?>> getResultClasses() {
    return Arrays.asList(String.class, JSONResult.class, Object.class);
  }

	@Override
	public <R> R getResults(Class<? extends R> type) {
		if(embedding == null)
			return null;
    if (type.equals(String.class)) {
      return (R)"Complete";
    } else if (type.equals(JSONResult.class)) {
        JSONResult res = () -> {
        StringBuilder builder = new StringBuilder();
        builder.append("[");
        for (String row: embedding.keySet()) {
          builder.append("{"+row+": ["+embedding.get(row)[0]+","+embedding.get(row)[1]+"]},");
        }
        return builder.substring(0, builder.length()-1)+"]";
      };
      return (R) res;
    }
		return (R) embedding;
	}
}
