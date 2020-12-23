package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

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
import edu.ucsf.rbvi.scNetViz.internal.view.ViewUtils;


// TODO: Consider random subsampling?
// TODO: Expose filter criteria
abstract public class AbstractEmbeddingTask extends AbstractTask implements ObservableTask{

	protected double[][] embedding;
	final ScNVManager manager;
		
	public AbstractEmbeddingTask(final ScNVManager manager) {
		this.manager = manager;
	}

	public void cancel() {
		cancelled = true;
	}

	public double[][] getResults() { return embedding; }

	public void scale(double[][] values) {
		double xMax = Double.MIN_VALUE;
		double xMin = Double.MAX_VALUE;
		double yMax = Double.MIN_VALUE;
		double yMin = Double.MAX_VALUE;
		for (int i = 0; i < values.length; i++) {
			xMin = Math.min(xMin, values[i][0]);
			xMax = Math.max(xMax, values[i][0]);
			yMin = Math.min(yMin, values[i][1]);
			yMax = Math.max(yMax, values[i][1]);
		}
		double xScale = xMax - xMin;
		double yScale = yMax - yMin;

		for (int i = 0; i < values.length; i++) {
			values[i][0] = (values[i][0]-xMin)/xScale;
			values[i][1] = (values[i][1]-yMin)/yScale;
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
        for (int row = 0; row < embedding.length; row++) {
          builder.append("{"+row+": ["+embedding[row][0]+","+embedding[row][1]+"]},");
        }
        return builder.substring(0, builder.length()-1)+"]";
      };
      return (R) res;
    }
		return (R) embedding;
	}
}
