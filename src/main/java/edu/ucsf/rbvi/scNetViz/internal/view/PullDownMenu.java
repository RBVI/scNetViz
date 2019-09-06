package edu.ucsf.rbvi.scNetViz.internal.view;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.Map;

import javax.swing.JComboBox;
import javax.swing.SwingUtilities;

import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskObserver;

import edu.ucsf.rbvi.scNetViz.internal.model.ScNVManager;

public class PullDownMenu extends JComboBox<String> {
	Map<String, Task> menuItems;
	final TaskObserver observer;
	final ScNVManager manager;
	final String title;
	boolean updating=false;

	public PullDownMenu(final ScNVManager manager, final String title,
	                    final Map<String, Task> menuItems, final TaskObserver observer) {
		super();
		this.manager = manager;
		this.menuItems = menuItems;
		this.observer = observer;
		this.title = title;
		init();
	}

	private void init() {
		addItem(title);
		for (String item: menuItems.keySet()) {
			addItem(item);
		}

		setFont(new Font("SansSerif", Font.PLAIN, 10));
		addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (updating) return;
				final String type = (String)getSelectedItem();
				if (menuItems.containsKey(type) && menuItems.get(type) != null) {
					final Task t = menuItems.get(type);
					SwingUtilities.invokeLater(new Runnable() {
								public void run() {
									manager.executeTasks(new TaskIterator(t), observer);
								}
					});
				}
				setSelectedIndex(0);
			}
		});
	}

	public void updateMenu(Map<String, Task> menuItems) {
		this.menuItems = menuItems;
		updating = true;
		removeAllItems();
		updating = false;
		init();
	}
}
