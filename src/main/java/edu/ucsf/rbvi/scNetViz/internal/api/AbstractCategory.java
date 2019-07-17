package edu.ucsf.rbvi.scNetViz.internal.api;

import java.io.File;
import java.io.InputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import org.apache.log4j.Logger;
// import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
// import org.apache.commons.math3.stat.ranking.NaNStrategy;
// import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.hipparchus.stat.inference.MannWhitneyUTest;
import org.hipparchus.stat.ranking.NaNStrategy;
import org.hipparchus.stat.ranking.TiesStrategy;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.work.TaskMonitor;

import edu.ucsf.rbvi.scNetViz.internal.api.Category;
import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.DoubleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.IntegerMatrix;
import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;
import edu.ucsf.rbvi.scNetViz.internal.api.StringMatrix;
import edu.ucsf.rbvi.scNetViz.internal.model.MatrixMarket;
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.SimpleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVWriter;
import edu.ucsf.rbvi.scNetViz.internal.utils.LogUtils;
import edu.ucsf.rbvi.scNetViz.internal.view.SortableTableModel;

public abstract class AbstractCategory extends SimpleMatrix implements Category {
	final protected Experiment experiment;

	Map<Object, Integer> sizes = null;
	Map<Object, double[]> means = null;
	Map<Object, int[]> countMap = null;
	Map<Object, double[]> mtdcMap = null;
	Map<Object, double[]> drMap = null;
	Map<Object, List<Integer>> catMap = null;

	protected final String name;

	protected int selectedRow = -1;
	protected int lastCategory = -1;
	protected int hdrCols = 1;
	protected BitSet excludeRows = null;
	protected int geneRows = 0;

	public AbstractCategory(final ScNVManager scManager, final Experiment exp,
	                        final String name, int nRows, int nCols) {
		super(scManager);
		super.nRows = nRows;
		super.nCols = nCols;
		this.experiment = exp;
		this.name = name;

		// XXX - should this be Matrix api?
		updateMatrixInfo();
	}

	@Override
	public int getSelectedRow() { return selectedRow; }

	@Override
	public void setSelectedRow(int selectedRow) { 
		if (lastCategory != selectedRow)
			lastCategory = -1;
		this.selectedRow = selectedRow; 
	}

	@Override
	public Map<Object, double[]> getMTDC(int category) {
		if (means != null && category == lastCategory)
			return mtdcMap;
		getMeans(category);
		return mtdcMap;
	}

	@Override
	public Map<Object, double[]> getMeans(int category) {
		updateMatrixInfo();
		if (means != null && category == lastCategory) {
			return means;
		}

		if (sizes == null || category != lastCategory) {
			getSizes(category);
		}

		means = new HashMap<>();
		countMap = new HashMap<>();
		mtdcMap = new HashMap<>();

		drMap = null;

		Matrix mtx = experiment.getMatrix();
		DoubleMatrix dMat = null;
		IntegerMatrix iMat = null;
		if (mtx instanceof DoubleMatrix) {
			dMat = (DoubleMatrix)mtx;
		} else if (mtx instanceof IntegerMatrix) {
			iMat = (IntegerMatrix)mtx;
		}


		int [] totalCount = new int[geneRows];
		Arrays.fill(totalCount, 0);

		for (Object key: catMap.keySet()) {
			// Skip over the unused categories -- they probably didn't pass QC
			if (key.toString().equals(Category.UNUSED_CAT))
				continue;

			calculateMeans(key, dMat, iMat, mtx.getNRows(), totalCount);
		}
		countMap.put("Total", totalCount); // Remember the total counts for each gene
		lastCategory = category;
		return means;
	}

	public void calculateMeans(Object key, DoubleMatrix dMat, IntegerMatrix iMat, int rowCount, int[] totalCount) {
		// System.out.println("Category: "+key);
		List<Integer> arrays = catMap.get(key);
		double[] catMean = new double[geneRows];
		int[] catCount = new int[geneRows];
		double[] catMTDC = new double[geneRows];

		int rowNumber = 0;
		for (int row = 0; row < rowCount; row++) {
			if (excludeRows.get(row))
				continue;

			double mean = 0.0;
			int foundCount = 0;
			for (Integer col: arrays) {
				double v = 0.0;
				if (dMat != null) {
					v = dMat.getDoubleValue(row, col);
				} else if (iMat != null) {
					v = (double)iMat.getIntegerValue(row, col);
				}
				// if (row == 0) System.out.println("v["+row+"]["+(col+1)+"] = "+v);
				if (!Double.isNaN(v)) {
					mean += v;
					foundCount++;
				}
			}
			catMean[rowNumber] = mean/(double)arrays.size();
			catCount[rowNumber] = foundCount;
			synchronized (totalCount) {
				totalCount[rowNumber] += foundCount;
			}
			catMTDC[rowNumber] = mean/(double)foundCount;
			rowNumber++;
		}
		means.put(key, catMean);
		countMap.put(key, catCount);
		mtdcMap.put(key, catMTDC);
	}
	

	@Override
	public Map<Object, int[]> getCounts(int category) {
		if (means != null && category == lastCategory)
			return countMap;
		getMeans(category);
		return countMap;
	}

	@Override
	public Map<Object, double[]> getDr(int category) {
		updateMatrixInfo();
		if (drMap != null && category == lastCategory)
			return drMap;
		if (means == null || category != lastCategory)
			getMeans(category);

		Matrix mtx = experiment.getMatrix();

		drMap = new HashMap<>();
		int totalAssays = mtx.getNCols();

		for (Object cat: means.keySet()) {
			calculateDr(cat, totalAssays);
		}
		return drMap;
	}

	public void calculateDr(Object cat, int totalAssays) {
		// System.out.println("Cat: "+cat);
		double dr[] = new double[geneRows];
		for (int row = 0; row < geneRows; row++) {
			double pct1 = (double)countMap.get(cat)[row]/(double)sizes.get(cat);
			double pct2 = (double)(countMap.get("Total")[row]-countMap.get(cat)[row])/(double)(totalAssays-sizes.get(cat));
			// dr[row] = Math.abs(pct1-pct2);
			dr[row] = Math.max(pct1,pct2);
		}
		drMap.put(cat, dr);
	}

	@Override
	public Map<Object, Integer> getSizes(int category) {
		if (sizes == null || category != lastCategory) {
			// This creates the sizes map as a by-product
			getUniqValues(category);
		}
		return sizes;
	}

	@Override
	public Map<Object, List<Integer>> getCatMap(int category) {
		if (sizes == null || category != lastCategory) {
			// This creates the sizes map as a by-product
			getUniqValues(category);
		}
		return catMap;
	}

	@Override
	public Map<Object, Map<String,double[]>> getLogGER(int category, double dDRthreshold, double log2FCCutoff) {
		if (means == null || category != lastCategory)
			getMeans(category);

		Map<Object, Map<String,double[]>> logGER = new HashMap<>();
		ExecutorService threadPool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors()-1);
		List<Callable <Map <Object, Map <String, double[]>>>> processes = new ArrayList<>();
		for (Object cat: means.keySet()) {
			processes.add(new GetLogGER(category, cat, dDRthreshold, log2FCCutoff));
		}

		try {
			List<Future<Map<Object, Map<String, double[]>>>> futures = threadPool.invokeAll(processes);
			for (Future<Map<Object, Map<String, double[]>>> future: futures) {
				Map<Object, Map<String, double[]>> ger = future.get();
				for (Object key: ger.keySet())
					logGER.put(key, ger.get(key));
			}
		} catch (Exception e) {}

		return logGER;

	}

	@Override
	public Map<String, double[]> getLogGER(int category, Object category1, double dDRthreshold, double log2FCCutoff) {
		if (means == null || category != lastCategory)
			getMeans(category);

		updateMatrixInfo();

		double[] logFC = new double[geneRows];
		double[] catMean = means.get(category1);
		double[] otherMean = new double[geneRows];
		double[] pValue = new double[geneRows];
		int[] thisCount = countMap.get(category1);
		int cat1Size = sizes.get(category1);
		int otherSize = 0;

		List<Integer> cat1Columns = catMap.get(category1);
		List<Integer> otherCatColumns = new ArrayList<>();

		MannWhitneyUTest mwTest = new MannWhitneyUTest(NaNStrategy.MINIMAL, TiesStrategy.AVERAGE);

		for (Object cat2: means.keySet()) {
			if (category1.equals(cat2))
				continue;

			otherCatColumns.addAll(catMap.get(cat2));

			int thisSize = sizes.get(cat2);
			double[] thisMean = means.get(cat2);
			int[] cat2Counts = countMap.get(cat2);
			otherSize += thisSize;

			// Calculate the mean of all other values
			for (int row = 0; row < geneRows; row++) {
				otherMean[row] += thisMean[row]*thisSize;
			}
		}

		// Now calculate our actual dr

		// System.out.println("cat1Columns.length = "+cat1Columns.size()+", otherCatColumns.length = "+otherCatColumns.size()+
		//                    " Matrix columns = "+experiment.getMatrix().getNCols());

		double log2 = Math.log(2.0);

		// Calculate the log2 fold change
		double[] dr = drMap.get(category1);
		for (int row = 0; row < geneRows; row++) {
			if (dr[row] < dDRthreshold) {
				logFC[row] = Double.NaN;
				pValue[row] = Double.NaN;
				continue;
			}

			logFC[row] = calcLogFC(catMean[row], cat1Size, otherMean[row]/otherSize, otherSize)/log2;

			// Calculate the MannWhitneyUTest
			if (Math.abs(logFC[row]) > log2FCCutoff)
				pValue[row] = mannWhitney(mwTest, row, cat1Columns, otherCatColumns);
			else
				pValue[row] = Double.NaN;
		}
		Map<String, double[]> returnMap = new HashMap<>();
		returnMap.put("logFC", logFC);
		returnMap.put("pValue", pValue);
		return returnMap;
	}

	@Override
	public double[] getLogGER(int category, Object category1, Object category2, 
	                          double dDRthreshold, double log2FCCutoff) {
		if (means == null || category != lastCategory)
			getMeans(category);

		updateMatrixInfo();

		double[] logFC = new double[geneRows];
		double[] catMean = means.get(category1);
		double[] otherMean = means.get(category2);
		double[] thisDr = drMap.get(category1);
		int cat1Size = sizes.get(category1);
		int cat2Size = sizes.get(category2);

		double log2 = Math.log(2.0);

		// Calculate the log2 fold change
		for (int row = 0; row < geneRows; row++) {
			if (thisDr[row] < dDRthreshold)
				continue;

			logFC[row] = calcLogFC(catMean[row], cat1Size, otherMean[row]/cat2Size, cat2Size)/log2;

			if (Math.abs(logFC[row]) <= log2FCCutoff)
				logFC[row] = Double.NaN;
		}
		return logFC;
	}

	public void saveFile(File file) throws IOException {
		CSVWriter.writeCSV(file, this, "\t");
	}

	double calcLogFC(double mean1, int size1, double mean2, int size2) {
		if (mean1 != 0.0 && mean2 != 0.0)
			return Math.log(mean1/mean2);

		if (mean2 == 0.0) {
			return Math.log(mean1/(1.0/size2));
		}
		return Math.log((1.0/size1)/mean2);
	}

	abstract public Object getValue(int row, int col);

	/**
	 * Note this returns an array of indices in the *experiment* matrix,
	 * not in the category matrix.  This makes it much easier later on
	 * to use this.
	 */
	protected int getUniqValues(int row) {
		// We'll use the column names to map to the experiment matrix
		Matrix mtx = experiment.getMatrix();
		HashMap<String, Integer> colLabelMap = new HashMap<>();
		List<String> colLabels = mtx.getColLabels();
		// System.out.println("nCols = "+nCols+", colLabels.size() = "+colLabels.size());
		int colIndex = 0;
		for (String str: colLabels) {
			colLabelMap.put(str, colIndex);
			colIndex++;
		}

		catMap = new HashMap<>();
		sizes = new HashMap<>();

		int tpmHeaderCols = 1;

		for (int col = hdrCols; col < nCols; col++) {
			// Get the corresponding mtx column first
			int mtxCol = mapColumn(col, tpmHeaderCols, colLabelMap);
			if (mtxCol < 0) continue;
			Object v = getValue(row, col);
			if (v != null && v.toString().length() == 0) {
				v = "Unspecified";
			}
			if (!catMap.containsKey(v)) {
				catMap.put(v, new ArrayList<>());
				sizes.put(v, 0);
			}
			catMap.get(v).add(mtxCol);
			sizes.put(v, sizes.get(v)+1);
		}

		if (colLabelMap.size() > 0) {
			// OK, now, any keys available in colLabelMap are unused.  Create a new category for them
			catMap.put(UNUSED_CAT, new ArrayList<>());
			sizes.put(UNUSED_CAT, 0);
			for (String key: colLabelMap.keySet()) {
				// System.out.println("Adding "+key+" to catMap");
				catMap.get(UNUSED_CAT).add(colLabelMap.get(key));
				sizes.put(UNUSED_CAT, sizes.get(UNUSED_CAT)+1);
			}
		}

		/*
		System.out.println("catMap: ");
		for (Object cat: catMap.keySet()) {
			System.out.println(cat+": "+catMap.get(cat));
		}
		*/
		return catMap.keySet().size();
	}

	private int mapColumn(int col, int tpmHeaders, Map<String, Integer> colLabelMap) {
		String lbl = getColumnLabel(col);
		if (colLabelMap.get(lbl) == null) {
			System.out.println("Can't find column label: "+lbl);
			return -1;
		}
		int mtxCol =  colLabelMap.get(lbl);
		colLabelMap.remove(lbl);
		return mtxCol;
	}

	protected double mannWhitney(MannWhitneyUTest test, int row, List<Integer> category1, 
	                             List<Integer> category2) {
		int assays = experiment.getMatrix().getNCols();

		double[] x = new double[category1.size()];
		double[] y = new double[category2.size()];

		Matrix mtx = experiment.getMatrix();
		DoubleMatrix dMat = null;
		IntegerMatrix iMat = null;
		if (mtx instanceof DoubleMatrix) {
			dMat = (DoubleMatrix)mtx;
		} else if (mtx instanceof IntegerMatrix) {
			iMat = (IntegerMatrix)mtx;
		}

		int index = 0;
		for (Integer col: category1) {
			double v = 0.0;
			if (dMat != null) {
				v = dMat.getDoubleValue(row, col);
			} else if (iMat != null) {
				v = (double)iMat.getIntegerValue(row, col);
			}
			x[index++] = v;
		}

		index = 0;
		for (Integer col: category2) {
			double v = 0.0;
			if (dMat != null) {
				v = dMat.getDoubleValue(row, col);
			} else if (iMat != null) {
				v = (double)iMat.getIntegerValue(row, col);
			}
			y[index++] = v;
		}

		double pValue = test.mannWhitneyUTest(x, y);
		return pValue;
	}

	protected int[][] getIntegerMatrix(boolean transpose) {
		if (transpose)
			return new int[nCols][nRows];
		else
			return new int[nRows][nCols];
	}

	protected String[][] getStringMatrix(boolean transpose) {
		if (transpose)
			return new String[nCols][nRows];
		else
			return new String[nRows][nCols];
	}

	protected void updateMatrixInfo() {
		if (experiment.getMatrix() == null) return;
		excludeRows = ((MatrixMarket)experiment.getMatrix()).findControls();
		geneRows = experiment.getMatrix().getNRows() - excludeRows.cardinality();
	}

	class GetLogGER implements Callable <Map <Object, Map <String, double[]>>> {
		int category;
		Object cat;
		double dDRthreshold;
		double log2FCCutoff;
		public GetLogGER(int category, Object cat, double dDRthreshold, double log2FCCutoff) {
			this.category = category;
			this.cat = cat;
			this.dDRthreshold = dDRthreshold;
			this.log2FCCutoff = log2FCCutoff;
		}

		public Map <Object, Map <String, double[]>> call() {
			Map<Object, Map<String, double[]>> returnMap = new HashMap<>();;
			returnMap.put(cat, getLogGER(category, cat, dDRthreshold, log2FCCutoff));
			return returnMap;
		}
	}

}
