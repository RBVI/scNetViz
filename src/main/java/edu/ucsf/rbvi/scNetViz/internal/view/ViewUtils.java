package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

public class ViewUtils {
	static int FONT_SIZE = 10;

	public static JButton addButton(ActionListener listener, String label, String command) {
		JButton button = new JButton(label);
		button.setFont(new Font("SansSerif", Font.PLAIN, FONT_SIZE));
		button.setAlignmentX(Component.CENTER_ALIGNMENT);
		button.setActionCommand(command);
		button.addActionListener(listener);
		return button;
	}

	public static JRadioButton addRadioButton(ActionListener listener, ButtonGroup group, 
	                                          String label, String command, boolean set) {
		JRadioButton button = new JRadioButton(label);
		button.setFont(new Font("SansSerif", Font.PLAIN, FONT_SIZE));
		button.setActionCommand(command);
		button.addActionListener(listener);
		button.setSelected(set);
		group.add(button);
		return button;
	}

	public static JTextField addLabeledField(JPanel container, String label, String defaultValue) {
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));
		panel.add(Box.createHorizontalGlue());

		JLabel jlabel = new JLabel(label);
		jlabel.setFont(new Font("SansSerif", Font.BOLD, FONT_SIZE));
		jlabel.setMaximumSize(new Dimension(100,35));
		panel.add(jlabel);

		JTextField field = new JTextField(defaultValue);
		field.setFont(new Font("SansSerif", Font.PLAIN, FONT_SIZE));
		field.setMinimumSize(new Dimension(50,25));
		field.setPreferredSize(new Dimension(50,25));
		field.setMaximumSize(new Dimension(50,35));
		panel.add(field);
		panel.add(Box.createRigidArea(new Dimension(15,0)));
		container.add(panel);
		return field;
	}

	public static JPanel addJCheckBox(ActionListener listener, String label, 
	                                  String command, boolean selected) {
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));
		panel.add(Box.createHorizontalGlue());

		JCheckBox checkBox = new JCheckBox(label, selected);
		checkBox.setHorizontalTextPosition(SwingConstants.LEFT);
		checkBox.setFont(new Font("SansSerif", Font.BOLD, FONT_SIZE));
		checkBox.setActionCommand(command);
		checkBox.addActionListener(listener);
		panel.add(checkBox);
		panel.add(Box.createRigidArea(new Dimension(15,0)));
		return panel;
	}

}
