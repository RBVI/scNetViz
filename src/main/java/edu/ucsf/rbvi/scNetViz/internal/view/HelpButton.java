package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.HashMap;
import java.util.Map;

import javax.swing.JButton;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class HelpButton extends JButton {
	final ScNVManager scManager;
	final String subPart;

	public HelpButton(final ScNVManager scManager, String subPart)
	{
		super("Help");
		this.scManager = scManager;
		this.subPart = subPart;

		setFont(new Font("SansSerif", Font.PLAIN, 10));
		addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Map<String, Object> args = new HashMap<>();
				args.put("id","scNetViz");
				args.put("title", "scNetViz Help");
				if (subPart != null) {
					args.put("url", "http://preview.rbvi.ucsf.edu/cytoscape/scNetViz/index.shtml#"+subPart);
				} else {
					args.put("url", "http://preview.rbvi.ucsf.edu/cytoscape/scNetViz/index.shtml");
				}
				scManager.executeCommand("cybrowser", "dialog", args, false);
			}
		});
	}
}
