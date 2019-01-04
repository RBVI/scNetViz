package edu.ucsf.rbvi.scNetViz.internal;

import static org.cytoscape.work.ServiceProperties.COMMAND;
import static org.cytoscape.work.ServiceProperties.COMMAND_DESCRIPTION;
import static org.cytoscape.work.ServiceProperties.COMMAND_EXAMPLE_JSON;
import static org.cytoscape.work.ServiceProperties.COMMAND_LONG_DESCRIPTION;
import static org.cytoscape.work.ServiceProperties.COMMAND_NAMESPACE;
import static org.cytoscape.work.ServiceProperties.COMMAND_SUPPORTS_JSON;
import static org.cytoscape.work.ServiceProperties.ID;
import static org.cytoscape.work.ServiceProperties.IN_MENU_BAR;
import static org.cytoscape.work.ServiceProperties.IN_TOOL_BAR;
import static org.cytoscape.work.ServiceProperties.INSERT_SEPARATOR_BEFORE;
import static org.cytoscape.work.ServiceProperties.LARGE_ICON_URL;
import static org.cytoscape.work.ServiceProperties.MENU_GRAVITY;
import static org.cytoscape.work.ServiceProperties.PREFERRED_MENU;
import static org.cytoscape.work.ServiceProperties.TITLE;
import static org.cytoscape.work.ServiceProperties.TOOL_BAR_GRAVITY;
import static org.cytoscape.work.ServiceProperties.TOOLTIP;

import java.util.Properties;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.application.swing.CytoPanelComponent2;
import org.cytoscape.io.BasicCyFileFilter;
import org.cytoscape.io.DataCategory;
import org.cytoscape.io.read.InputStreamTaskFactory;
import org.cytoscape.io.util.StreamUtil;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.work.TaskFactory;

import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.service.util.AbstractCyActivator;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.osgi.framework.BundleContext;
import org.osgi.framework.ServiceReference;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;
import edu.ucsf.rbvi.scNetViz.internal.sources.file.FileSource;
import edu.ucsf.rbvi.scNetViz.internal.tasks.ListExperimentsTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.tasks.SelectTaskFactory;
import edu.ucsf.rbvi.scNetViz.internal.tasks.SettingsTaskFactory;

public class CyActivator extends AbstractCyActivator {

	public CyActivator() {
		super();
	}

	public void start(BundleContext bc) {
		final StreamUtil streamUtil = getService(bc, StreamUtil.class);
		final CyServiceRegistrar serviceRegistrar = getService(bc, CyServiceRegistrar.class);

		final ScNVManager scNVManager = new ScNVManager(serviceRegistrar);

		// Commands
		//
		// Source-specific commands (registered by source)
		// *scnetviz load experiment file file=yyyy
		// *scnetviz load category file file=yyyy accession=yyyy experiment=exp
		// *scnetviz load experiment gxa accession=yyyy
		// *scnetviz list gxa entries
		//
		// Manager commands
		// *scnetviz list experiments
		// scnetviz remove experiment accession=yyyy
		// scnetviz show experiment accession=yyyy
		// scnetviz export experiment accession=yyyy
		//
		// Category commands
		// scnetviz export category accession=yyyy category=zzzz
		// scnetviz calculate diffexp accession=yyyy category=zzzz categoryrow=nnn logfc=ddd min.pct=ppp
		//
		// Differential expression commands
		// scnetviz export diffexp accesssion=yyyy
		// scnetviz create networks accession=yyyy category=zzzz categoryrow=nnn logfc=ddd min.pct=ppp pvalue=ddd log2fc=nnn topgenes=nnn maxgenes=nnn
		// scnetviz show experiment accession=yyyy
		// scnetviz show violin accession=yyyy
		// scnetviz show heatmap accession=yyyy

		// Register our sources
		scNVManager.addSource(new GXASource(scNVManager));
		scNVManager.addSource(new FileSource(scNVManager));

		{
			ListExperimentsTaskFactory list = new ListExperimentsTaskFactory(scNVManager);
			Properties props = new Properties();
			props.setProperty(COMMAND_NAMESPACE, "scnetviz");
			props.setProperty(COMMAND, "list experiments");
			props.setProperty(COMMAND_DESCRIPTION, "List the currently loaded experiments");
			props.setProperty(COMMAND_LONG_DESCRIPTION, "");
			props.setProperty(COMMAND_SUPPORTS_JSON, "true");
			props.setProperty(COMMAND_EXAMPLE_JSON, "{}");
			scNVManager.registerService(list, TaskFactory.class, props);
		}

		{
			SelectTaskFactory select = new SelectTaskFactory(scNVManager);
			Properties props = new Properties();
			props.setProperty(COMMAND_NAMESPACE, "scnetviz");
			props.setProperty(COMMAND, "select");
			props.setProperty(COMMAND_DESCRIPTION, "Select genes or assays in current tables");
			props.setProperty(COMMAND_LONG_DESCRIPTION, "");
			props.setProperty(COMMAND_SUPPORTS_JSON, "true");
			props.setProperty(COMMAND_EXAMPLE_JSON, "{}");
			scNVManager.registerService(select, TaskFactory.class, props);
		}


		{
			SettingsTaskFactory settings = new SettingsTaskFactory(scNVManager);
			Properties props = new Properties();
			props.put(TITLE, "Settings");
			props.put(PREFERRED_MENU, "Apps.scNetViz");
			props.setProperty(IN_TOOL_BAR, "FALSE");
			props.setProperty(IN_MENU_BAR, "TRUE");
			scNVManager.registerService(settings, TaskFactory.class, props);
		}

		/*
		{
			// This is for the basic reader.  Note that we'll also load a more advanced one below
			final BasicCyFileFilter mtxFileFilter = new BasicCyFileFilter(new String[] { "mtx" },
			                              new String[] { "application/mtx" }, "MTX", DataCategory.TABLE, streamUtil);
			final MTXReaderTaskFactory mtxReaderFactory = new MTXReaderTaskFactory(mtxFileFilter, scNVManager);
	
			Properties mtxReaderProps = new Properties();
			mtxReaderProps.put(ID, "mtxTableReaderFactory");
			registerService(bc, mtxReaderFactory, InputStreamTaskFactory.class, mtxReaderProps);
	
			Properties mtxImporterProps = new Properties();
			mtxImporterProps.setProperty(PREFERRED_MENU, "Apps.MTXImporter");
			mtxImporterProps.setProperty(TITLE, "Import MTX files");
			registerService(bc, mtxReaderFactory, TaskFactory.class, mtxImporterProps);
		}
		*/

	}
}
