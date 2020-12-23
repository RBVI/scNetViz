package edu.ucsf.rbvi.scNetViz.internal.tasks;

import java.io.File;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.ContainsTunables;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.ProvidesTitle;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;
import org.cytoscape.work.Tunable;
import org.cytoscape.work.util.ListSingleSelection;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.model.DifferentialExpression;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.view.ViewUtils;


// TODO: Consider random subsampling?
// TODO: Expose filter criteria
public class ShowDiffPlotTask extends AbstractTask {
  final ScNVManager manager;
  DifferentialExpression diffExp;

	@Tunable (description="Experiment to show plot for", context="nogui")
	public ListSingleSelection<String> accession = null;

	@Tunable (description="Type of plot to show", context="nogui", required=true)
	public ListSingleSelection<String> type = new ListSingleSelection<>("Heatmap","Violin");

  @Tunable (description="Gene",context="nogui",
            longDescription="Optional gene to use to narrow violin plot")
  public String gene;

  @Tunable (description="Positive only",context="nogui",
            longDescription="Only show genes with increased expression")
  public boolean posOnly = false;

  @Tunable (description="FDR cutoff",context="nogui",
            longDescription="FDR cutoff for heatmap")
  public double fdrCutoff = 0.5d;

  @Tunable (description="Log2FC cutoff",context="nogui",
            longDescription="Log 2 fold change (log2FC) for heatmap")
  public double log2FCCutoff = 1.0d;

	public ShowDiffPlotTask(final ScNVManager manager) {
		List<String> accessions = new ArrayList<String>(manager.getExperimentAccessions());
		accession = new ListSingleSelection<>(new ArrayList<String>(manager.getExperimentAccessions()));
    this.manager = manager;
	}

	public void run(TaskMonitor monitor) {
		// Get the experiment
		Experiment exp = manager.getExperiment(accession.getSelectedValue());

    diffExp = exp.getDiffExp();
    if (diffExp == null) {
      monitor.showMessage(TaskMonitor.Level.ERROR, "No differential expression calculated");
      return;
    }

    if (type.getSelectedValue().equals("Heatmap"))
      showHeatMap(exp);
    else if (type.getSelectedValue().equals("Violin"))
      showViolin(exp, gene);
	}

  void showHeatMap(Experiment exp) {
		Task heatTask = new HeatMapTask(manager, diffExp.getCategory(), diffExp, fdrCutoff, 
                                    log2FCCutoff, posOnly, -1, null);
		insertTasksAfterCurrentTask(heatTask);
  }

  void showViolin(Experiment exp, String gene) {
    if (gene == null) {
      Task vdeTask = new ViolinDiffExpTask(manager, diffExp.getCategory(), diffExp, posOnly);
			insertTasksAfterCurrentTask(vdeTask);
    } else {
      int row = 0;
      int rowNumber = -1;
      for (String rowName: exp.getMatrix().getRowLabels()) {
        if (rowName.equalsIgnoreCase(gene)) {
          rowNumber = row;
          break;
        }
        row++;
      }
			int selectedCategory = diffExp.getCategoryRow();
			Task vgeneTask = new ViolinGeneTask(manager, diffExp.getCategory(), selectedCategory, rowNumber);
			insertTasksAfterCurrentTask(vgeneTask);
    }
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
      for (String rowName: exp.getMatrix().getRowLabels()) {
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
