package edu.ucsf.rbvi.scNetViz.internal.sources.gxa.view;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;

import java.io.File;

import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import org.cytoscape.application.swing.CytoPanelComponent2;
import org.cytoscape.application.swing.CytoPanelName;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXAMetadata;
import edu.ucsf.rbvi.scNetViz.internal.sources.gxa.GXASource;

public class GXAEntryFrame extends JFrame {
	final ScNVManager scNVManager;
	JScrollPane scrollPane = null;
	GXAEntryTable gxaEntryTable;
	GXASource gxaSource;
	JCheckBox loadOnly;
	Icon ico = null;

	public GXAEntryFrame (final ScNVManager scNVManager, final GXASource gxaSource) {
		super();
		this.scNVManager = scNVManager;
		this.gxaSource = gxaSource;
		this.setLayout(new BorderLayout());
		this.setTitle("Single Cell Expression Atlas Browser");
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		init();
	}

	public void init() {
		if (scrollPane == null) {
			if (gxaSource.getMetadata().size() == 0) {
				gxaSource.loadGXAEntries(null);
			}

			JPanel bottomPanel = new JPanel(new BorderLayout());
			bottomPanel.add(new JLabel("<html><i>Double click on experiment row to launch network analysis</i></html>", SwingConstants.LEFT), 
			             BorderLayout.WEST);
			loadOnly = new JCheckBox("<html><span style=\"font-size: 80%\">Load experiment only (don't analyze)</span></html>");
			loadOnly.setHorizontalTextPosition(SwingConstants.LEFT);
			loadOnly.setSelected(false);
			bottomPanel.add(loadOnly, BorderLayout.EAST);

			this.add(bottomPanel, BorderLayout.SOUTH);

			GXAEntryTableModel tableModel = new GXAEntryTableModel(gxaSource.getMetadata());
			gxaEntryTable = new GXAEntryTable(scNVManager, gxaSource, tableModel);
			scrollPane = new JScrollPane(gxaEntryTable);
			scrollPane.setPreferredSize(new Dimension(900, 500));
			this.add(scrollPane, BorderLayout.CENTER);
			this.pack();
			this.setVisible(true);
		}
	}

	public boolean isLoadOnly() {
		return loadOnly.isSelected();
	}

	@Override
	public String getTitle() {
		return "Single Cell Expression Atlas Browser";
	}

}
