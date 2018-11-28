package edu.ucsf.rbvi.scNetViz.internal.api;

import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.model.SimpleMatrix;
import edu.ucsf.rbvi.scNetViz.internal.utils.CSVReader;
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
	protected int hdrCols = 1;

	public AbstractCategory(final ScNVManager scManager, final Experiment exp,
	                        final String name, int nRows, int nCols) {
		super(scManager);
		super.nRows = nRows;
		super.nCols = nCols;
		this.experiment = exp;
		this.name = name;
	}

	@Override
	public Map<Object, double[]> getMTDC(int category) {
		if (means != null && category == selectedRow)
			return mtdcMap;
		getMeans(category);
		return mtdcMap;
	}

	@Override
	public Map<Object, double[]> getMeans(int category) {
		if (means != null && category == selectedRow)
			return means;

		if (sizes == null || category != selectedRow) {
			getSizes(category);
		}

		means = new HashMap<>();
		countMap = new HashMap<>();
		mtdcMap = new HashMap<>();

		Matrix mtx = experiment.getMatrix();
		DoubleMatrix dMat = null;
		IntegerMatrix iMat = null;
		if (mtx instanceof DoubleMatrix) {
			dMat = (DoubleMatrix)mtx;
		} else if (mtx instanceof IntegerMatrix) {
			iMat = (IntegerMatrix)mtx;
		}

		int [] totalCount = new int[mtx.getNRows()];
		Arrays.fill(totalCount, 0);

		for (Object key: catMap.keySet()) {
			List<Integer> arrays = catMap.get(key);
			double[] catMean = new double[mtx.getNRows()];
			int[] catCount = new int[mtx.getNRows()];
			double[] catMTDC = new double[mtx.getNRows()];
			for (int row = 0; row < mtx.getNRows(); row++) {
				double mean = 0.0;
				int foundCount = 0;
				for (Integer col: arrays) {
					double v = 0.0;
					if (dMat != null) {
						v = dMat.getDoubleValue(row, col);
					} else if (iMat != null) {
						v = (double)iMat.getIntegerValue(row, col);
					}
					if (!Double.isNaN(v)) {
						mean += v;
						foundCount++;
					}
				}
				catMean[row] = mean/(double)arrays.size();
				catCount[row] = foundCount;
				totalCount[row] += foundCount;
				catMTDC[row] = mean/foundCount;
			}
			means.put(key, catMean);
			countMap.put(key, catCount);
			mtdcMap.put(key, catMTDC);
		}
		countMap.put("Total", totalCount); // Remember the total counts for each gene
		return means;
	}

	@Override
	public Map<Object, int[]> getCounts(int category) {
		if (means != null && category == selectedRow)
			return countMap;
		getMeans(category);
		return countMap;
	}

	@Override
	public Map<Object, double[]> getDr(int category) {
		if (drMap != null && category == selectedRow)
			return drMap;
		if (means == null || category != selectedRow)
			getMeans(category);

		drMap = new HashMap<>();
		for (Object cat: means.keySet()) {
			double dr[] = new double[experiment.getMatrix().getNRows()];
			for (int row = 0; row < experiment.getMatrix().getNRows(); row++) {
				double pct1 = (double)countMap.get(cat)[row]/(double)sizes.get(cat);
				double pct2 = (double)(countMap.get("Total")[row]-countMap.get(cat)[row])/(double)(nCols-sizes.get(cat));
				dr[row] = Math.abs(pct1-pct2);
			}
			drMap.put(cat, dr);
		}
		return drMap;
	}

	@Override
	public Map<Object, Integer> getSizes(int category) {
		if (sizes == null || category != selectedRow) {
			// This creates the sizes map as a by-product
			getUniqValues(category);
		}
		return sizes;
	}

	@Override
	public Map<Object, Map<String,double[]>> getLogGER(int category, double dDRthreshold) {
		if (means == null || category == selectedRow)
			getMeans(category);

		Map<Object, Map<String,double[]>> logGER = new HashMap<>();
		for (Object cat: means.keySet()) {
			logGER.put(cat, getLogGER(category, cat, dDRthreshold));
		}

		return logGER;

	}

	@Override
	public Map<String, double[]> getLogGER(int category, Object category1, double dDRthreshold) {
		if (means == null || category == selectedRow)
			getMeans(category);

		int geneRows = experiment.getMatrix().getNRows();
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

			logFC[row] = Math.log(catMean[row]/(otherMean[row]/otherSize))/log2;

			// Calculate the MannWhitneyUTest
			pValue[row] = mannWhitney(mwTest, row, cat1Columns, otherCatColumns);
		}
		Map<String, double[]> returnMap = new HashMap<>();
		returnMap.put("logFC", logFC);
		returnMap.put("pValue", pValue);
		return returnMap;
	}

	@Override
	public double[] getLogGER(int category, Object category1, Object category2, double dDRthreshold) {
		if (means == null || category == selectedRow)
			getMeans(category);

		int geneRows = experiment.getMatrix().getNRows();
		double[] logFC = new double[geneRows];
		double[] catMean = means.get(category1);
		double[] otherMean = means.get(category2);
		double[] thisDr = drMap.get(category1);
		int otherSize = 0;

		double log2 = Math.log(2.0);

		// Calculate the log2 fold change
		for (int row = 0; row < geneRows; row++) {
			if (thisDr[row] < dDRthreshold)
				continue;
			logFC[row] = Math.log(catMean[row]/(otherMean[row]))/log2;
		}
		return logFC;
	}

	abstract public Object getValue(int row, int col);

	protected int getUniqValues(int row) {
		catMap = new HashMap<>();
		sizes = new HashMap<>();
		for (int col = 0; col < nCols; col++) {
			Object v = getValue(row, col);
			if (!catMap.containsKey(v)) {
				catMap.put(v, new ArrayList<>());
				sizes.put(v, -1);
			}
			catMap.get(v).add(col);
			sizes.put(v, sizes.get(v)+1);
		}
		return catMap.keySet().size();
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

}
