package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cytoscape.service.util.CyServiceRegistrar;

import edu.ucsf.rbvi.scNetViz.internal.api.Matrix;

public class MatrixManager {
	final Map<String, Matrix> matrixMap;
	final ScNVManager manager;

	public MatrixManager(ScNVManager manager) {
		matrixMap = new HashMap<>();
		this.manager = manager;
	}

	public void addMatrix(String name, Matrix matrix) {
		matrixMap.put(name, matrix);
	}

	public Matrix getMatrix(String name) {
		if (matrixMap.containsKey(name)) return matrixMap.get(name);
		return null;
	}

	public Set<String> getMatrixNames() {
		return matrixMap.keySet();
	}

	public <S> S getService(Class<S> serviceClass) {
		return manager.getService(serviceClass);
	}
}
