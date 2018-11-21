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
	Map<Object, double[]> drMap = null;
	Map<Object, double[]> mtdcMap = null;
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
		drMap = new HashMap<>();
		mtdcMap = new HashMap<>();

		Matrix mtx = experiment.getMatrix();
		DoubleMatrix dMat = null;
		IntegerMatrix iMat = null;
		if (mtx instanceof DoubleMatrix) {
			dMat = (DoubleMatrix)mtx;
		} else if (mtx instanceof IntegerMatrix) {
			iMat = (IntegerMatrix)mtx;
		}

		for (Object key: catMap.keySet()) {
			List<Integer> arrays = catMap.get(key);
			double[] catMean = new double[mtx.getNRows()];
			double[] catDr = new double[mtx.getNRows()];
			double[] catMTDC = new double[mtx.getNRows()];
			for (int row = 0; row < mtx.getNRows(); row++) {
				double mean = 0.0;
				double dr = 0.0;
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
						dr ++;
						foundCount++;
					}
				}
				catMean[row] = mean/(double)arrays.size();
				catDr[row] = dr/(double)arrays.size();
				catMTDC[row] = mean/foundCount;
			}
			means.put(key, catMean);
			drMap.put(key, catDr);
			mtdcMap.put(key, catMTDC);
		}
		return means;
	}

	@Override
	public Map<Object, double[]> getDr(int category) {
		if (means != null && category == selectedRow)
			return drMap;
		getMeans(category);
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

	abstract public Object getValue(int row, int col);

	protected int getUniqValues(int row) {
		catMap = new HashMap<>();
		sizes = new HashMap<>();
		for (int col = hdrCols; col < nCols; col++) {
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

}
