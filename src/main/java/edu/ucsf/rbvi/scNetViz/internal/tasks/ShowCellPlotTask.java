package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ContainsTunables;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.utils.HTTPUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.ViewUtils;


// TODO: Consider random subsampling?
// TODO: Expose filter criteria
public class ShowCellPlotTask extends AbstractTask {
  final ScNVManager manager;

	@Tunable (description="Experiment to show plot for")
	public ListSingleSelection<String> accession = null;

  @Tunable (description="Category",context="nogui",
            longDescription="Category to use for coloring the plot")
  public String category;

  @Tunable (description="Category row",context="nogui",
            longDescription="Category row to use for coloring the plot")
  public int categoryRow;

  @Tunable (description="Gene",context="nogui",
            longDescription="Gene to use for coloring the plot")
  public String gene;

	public ShowCellPlotTask(final ScNVManager manager) {
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
    this.manager = manager;
	}

	public void run(TaskMonitor monitor) {
		// Get the experiment
		Experiment exp = manager.getExperiment(accession.getSelectedValue());
		// Get the MatrixMarket matrix
		Matrix mtx = exp.getMatrix();
		if (!(mtx instanceof MatrixMarket)) {
			monitor.showMessage(TaskMonitor.Level.ERROR,"Matrix must be of type MatrixMarket");
			return;
		}
		MatrixMarket mmtx = (MatrixMarket)mtx;

    if (exp.getPlotType() == null) {
			monitor.showMessage(TaskMonitor.Level.ERROR,"Plot must be calculated first");
			return;
    }

    showPlot(exp, category, categoryRow, gene);
	}

  protected void showPlot(Experiment exp, String categoryName, int categoryRow, String gene) {
    Category cat = null;
    if (categoryName != null) {
		  List<Category> categories = exp.getCategories();
		  if (categories != null) {
        for (Category c : categories) {
          if (c.toString().equalsIgnoreCase(categoryName)) {
            cat = c;
            break;
          }
        }
      }
		}

    int rowNumber = -1;
    if (gene != null) {
      int row = 0;
      for (String rowName: exp.getMatrix().getRowLabels(0)) {
        if (rowName.equalsIgnoreCase(gene)) {
          rowNumber = row;
          break;
        }
        row++;
      }
    }
    ViewUtils.showPlot(manager, exp, cat, categoryRow, rowNumber);
  }
}
