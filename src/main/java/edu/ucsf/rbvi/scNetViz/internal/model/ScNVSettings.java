package edu.ucsf.rbvi.scNetViz.internal.model;

import java.util.Properties;

import org.cytoscape.property.AbstractConfigDirPropsReader;
import org.cytoscape.property.CyProperty;

public class ScNVSettings extends AbstractConfigDirPropsReader {
	public enum SETTING {
		MAX_GENES("maxGenes", "200"),
		DE_FC_CUTOFF("deFcCutoff", "0.5"),
		DE_MIN_PCT_CUTOFF("deMinPctCutoff", "10"),
		NET_PV_CUTOFF("netPvCutoff","0.05"),
		NET_FC_CUTOFF("netFcCutoff","1.00"),
		POSITIVE_ONLY("positiveOnly","false"),
		DONT_ANALYZE("dontAnalyze","true"),
		HEATMAP_COUNT("heatMapCount","20");

		String name;
		String value;
		SETTING(String name, String value) {
			this.name = name;
			this.value = value;
		}

		public String getName() { return name; }
		public String getValue() { return value; }
	}

	public ScNVSettings() {
		super("scNetViz", "scNetViz.props", CyProperty.SavePolicy.SESSION_FILE_AND_CONFIG_DIR);

		// Initialize if our properties are not set
		for (SETTING setting: SETTING.values()) {
			if (!getProperties().containsKey(setting.getName()))
				getProperties().setProperty(setting.getName(), setting.getValue());
		}
	}

	public String getSetting(SETTING setting) {
		return getProperties().getProperty(setting.getName());
	}

	public void setSetting(SETTING setting, String value) {
		getProperties().setProperty(setting.getName(), value);
	}
}
