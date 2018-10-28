package edu.ucsf.rbvi.scNetViz.internal.view;

import javax.swing.JLabel;

import edu.ucsf.rbvi.scNetViz.internal.api.Experiment;
import edu.ucsf.rbvi.scNetViz.internal.api.Metadata;
import edu.ucsf.rbvi.scNetViz.internal.api.Source;

public class ExperimentLabel extends JLabel {

	public ExperimentLabel(final Experiment experiment) {
		super(createLabelText(experiment.getSource(), experiment.getMetadata()));
	}

	private static String createLabelText(Source source, Metadata meta) {
		String title = "<html><span style=\"font-size:10pt\">"+source.toString()+" <b>"+meta.get(Metadata.ACCESSION)+"</b></span><br/>";
		title += "<span style=\"font-size:10pt\"><i>"+meta.get(Metadata.DESCRIPTION)+"</i></span></html>";
		return title;
	}
}
