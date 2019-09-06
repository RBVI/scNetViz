package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Color;
import java.awt.Dimension;
import javax.swing.JTextPane;
import javax.swing.SwingConstants;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;

public class ExperimentLabel extends JTextPane {
	int width;

	public ExperimentLabel(final Experiment experiment, final Color background) {
		this(experiment, background, 250);
	}

	public ExperimentLabel(final Experiment experiment, final Color background, int width) {
		super();
		this.width = width;
		setContentType("text/html");
		setEditable(false);
		setBackground(background);
		setBorder(null);
		// putClientProperty(JEditorPane.HONOR_DISPLAY_PROPERTIES, true);
		setMaximumSize(new Dimension(250,75));

		setText(createLabelText(experiment.getSource(), experiment.getMetadata(),background, width));
	}

	public void updateText(Experiment experiment) {
		setText(createLabelText(experiment.getSource(), experiment.getMetadata(),getBackground(), width));
	}

	private static String createLabelText(Source source, Metadata meta, Color background, int width) {
		String bg = "background-color:rgb("+background.getRed()+
		            ","+background.getGreen()+","+background.getBlue()+")";
		String title = "<html><body style=\"text-align:center;font-size:10pt;"+bg+"\">";
		title += "<div style=\"width:"+width+"px\">";
		title += source.toString()+" <b>"+meta.get(Metadata.ACCESSION)+"</b><br/>";
		title += "<i>"+meta.get(Metadata.DESCRIPTION)+"</i>";
		title += "</div></body></html>";
		return title;
	}
}
