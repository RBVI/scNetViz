package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Color;
import java.awt.Dimension;
import javax.swing.JTextPane;
import javax.swing.SwingConstants;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;

public class ExperimentLabel extends JTextPane {

	public ExperimentLabel(final Experiment experiment, final Color background) {
		super();
		setContentType("text/html");
		setEditable(false);
		setBackground(background);
		setBorder(null);
		// putClientProperty(JEditorPane.HONOR_DISPLAY_PROPERTIES, true);
		setMaximumSize(new Dimension(250,75));

		setText(createLabelText(experiment.getSource(), experiment.getMetadata(),background));
	}

	public void updateText(Experiment experiment) {
		setText(createLabelText(experiment.getSource(), experiment.getMetadata(),getBackground()));
	}

	private static String createLabelText(Source source, Metadata meta, Color background) {
		String bg = "background-color:rgb("+background.getRed()+
		            ","+background.getGreen()+","+background.getBlue()+")";
		String title = "<html><body style=\"text-align:center;font-size:10pt;"+bg+"\">";
		title += "<div style=\"width:250px\">";
		title += source.toString()+" <b>"+meta.get(Metadata.ACCESSION)+"</b><br/>";
		title += "<i>"+meta.get(Metadata.DESCRIPTION)+"</i>";
		title += "</div></body></html>";
		return title;
	}
}
