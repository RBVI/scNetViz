package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cytoscape.service.util.CyServiceRegistrar;

public class MTXManager {
	final Map<String, MatrixMarket> mtxMap;
	final CyServiceRegistrar registrar;

	public MTXManager(CyServiceRegistrar registrar) {
		mtxMap = new HashMap<>();
		this.registrar = registrar;
	}

	public void addMatrix(String name, MatrixMarket mtx) {
		mtxMap.put(name, mtx);
	}

	public MatrixMarket getMatrix(String name) {
		if (mtxMap.containsKey(name)) return mtxMap.get(name);
		return null;
	}

	public Set<String> getMatrixNames() {
		return mtxMap.keySet();
	}

	public CyServiceRegistrar getRegistrar() { return registrar; }
	public <S> S getService(Class<S> serviceClass) { return registrar.getService(serviceClass); }
}
