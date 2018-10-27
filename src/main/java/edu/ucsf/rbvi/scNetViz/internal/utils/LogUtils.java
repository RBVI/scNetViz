package edu.ucsf.rbvi.scNetViz.internal.utils;

import org.apache.log4j.Logger;

import org.cytoscape.application.CyUserLog;
import org.cytoscape.work.TaskMonitor;

public class LogUtils {
	static Logger logger = Logger.getLogger(CyUserLog.NAME);

	public static void info(String message) {
		log(null, TaskMonitor.Level.INFO, message);
	}

	public static void warn(String message) {
		log(null, TaskMonitor.Level.WARN, message);
	}

	public static void error(String message) {
		log(null, TaskMonitor.Level.ERROR, message);
	}

	public static void log(TaskMonitor taskMonitor, TaskMonitor.Level level, String message) {
		if (taskMonitor != null) {
			taskMonitor.showMessage(level, message);
			return;
		}
		switch (level) {
			case ERROR:
				logger.error(message);
				break;
			case INFO:
				logger.info(message);
				break;
			case WARN:
				logger.warn(message);
				break;
		}
	}
}
