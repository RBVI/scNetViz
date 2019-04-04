package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Dimension;
import javax.swing.JTextPane;
import javax.swing.SwingConstants;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;

public class ExperimentLabel extends JTextPane {

	public ExperimentLabel(final Experiment experiment) {
		super();
		setContentType("text/html");
		setEditable(false);
		setBackground(null);
		setBorder(null);
		// putClientProperty(JEditorPane.HONOR_DISPLAY_PROPERTIES, true);
		setMaximumSize(new Dimension(250,75));
						
		setText(createLabelText(experiment.getSource(), experiment.getMetadata()));
	}

	public void updateText(Experiment experiment) {
		setText(createLabelText(experiment.getSource(), experiment.getMetadata()));
	}

	private static String createLabelText(Source source, Metadata meta) {
		String title = "<html><body style=\"text-align: center;font-size:10pt\">";
		title += "<div style=\"width:250px\">";
		title += source.toString()+" <b>"+meta.get(Metadata.ACCESSION)+"</b><br/>";
		title += "<i>"+meta.get(Metadata.DESCRIPTION)+"</i>";
		title += "</div></body></html>";
		/*
		String title = "<html><body style=\"text-align:center;font-size:10pt\">"+
						source.toString()+" <b>"+meta.get(Metadata.ACCESSION)+"</b>";
		title += "<div style=\"text-align:center;position:absolute;left:50%;margin: -125px 0 0 20px;width:250px;\"><i>"+meta.get(Metadata.DESCRIPTION)+"</i></div></body></html>";
		*/
		return title;
	}
}
