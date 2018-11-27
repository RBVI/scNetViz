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
		String title = "<html><span style=\"font-size:10pt\">"+
						source.toString()+" <b>"+meta.get(Metadata.ACCESSION)+"</b></span>";
		title += "<div style=\"font-size:10pt;width:250px\"><i>"+meta.get(Metadata.DESCRIPTION)+"</i></div></html>";
		return title;
	}
}
