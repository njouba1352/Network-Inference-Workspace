// TreeDataViewer.java
// the 'condition chooser' gui wraps the DataMatrixViewer
//-----------------------------------------------------------------------------------------------------
// RCSid = '$Revision: 3545 $   $Date: 2005/04/13 16:52:58 $'
//-----------------------------------------------------------------------------------------------------
/*
 * Copyright (C) 2006 by Institute for Systems Biology,
 * Seattle, Washington, USA.  All rights reserved.
 *
 * This source code is distributed under the GNU Lesser
 * General Public License, the text of which is available at:
 *   http://www.gnu.org/copyleft/lesser.html
 */

package org.systemsbiology.gaggle.experiment.gui;
//-----------------------------------------------------------------------------------------------------

import java.util.*;
import java.util.List;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import javax.swing.tree.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileFilter;

import org.systemsbiology.gaggle.experiment.metadata.*;
import org.systemsbiology.gaggle.experiment.gui.movie.*;
import org.systemsbiology.gaggle.experiment.datamatrix.*;
import org.systemsbiology.gaggle.experiment.readers.DataMatrixFileReader;

import org.systemsbiology.gaggle.util.MiscUtil;

//-----------------------------------------------------------------------------------------------------
public class TreeDataViewer implements TreeSelectionListener, WindowListener, java.io.Serializable {


    List<String> perst = new ArrayList<String>();

    // ArrayList dataDirectories;
    ArrayList<TreePath> currentTreeSelection = new ArrayList<TreePath>();
    DataMatrixViewer dataMatrixViewer;
    MetaDataNavigator experimentNavigator;
    ArrayList<MetaData> currentlySelectedExperiments = new ArrayList<MetaData>();
    ArrayList experiments = new ArrayList();
    JFrame jframe;
    JMenuBar menubar;
    JToolBar toolbar, extraToolbar;
    JPanel toolbarEtcPanel;
    JScrollPane treePane;
    JPanel resultsPanel;

    JSplitPane outerSplitPane;
    JTree treeWidget;
    JButton newRepositoryButton;
    JButton saveAllButton;
    JButton loadSelectedConditionsButton;
    JButton showHideWorkFramesButton;

    JButton followExperimentLinksButton;
    HashMap selectedExperimentInfoUrls;   // a hash of linkType -> linkUrl

    HashSet<ExitListener> exitListeners = new HashSet<ExitListener>();

    boolean workFramesVisible = true;
    JTextField conditionCounterTextField;
    boolean dataBrowserHasEntireWindow;

    Settings settings = Settings.getInstance();
    int selectedConditionCount;
    MovieController movieController;
    String repository;
    private ExperimentDetailsAction experimentDetailsAction;
    private CloseAllTabsAction closeAllTabsAction;

    boolean treeCollapsed = true;
    String preSelectedPaths;

    //-----------------------------------------------------------------------------------------------------
    public TreeDataViewer() {
        this(null, null);
    }

    public TreeDataViewer(String repository) {
        this(repository, null);
    } // ctor


    public TreeDataViewer(String repository, String preSelectedPaths) {
        //System.out.println ("TDV ctor: " + repository);
        this.repository = repository;
        this.preSelectedPaths = preSelectedPaths;
        //currentTreeSelection = new ArrayList ();
        //currentlySelectedExperiments = new ArrayList ();
        experiments = new ArrayList();
        createFrame();
        initializeSelectionWidgetsAndVariables();
        try {
            dataMatrixViewer = new DataMatrixViewer(experimentNavigator);
            resultsPanel.add(dataMatrixViewer, BorderLayout.CENTER);
            resultsPanel.paintAll(resultsPanel.getGraphics());
        }
        catch (Exception ex0) {
            ex0.printStackTrace();
        }
        if (preSelectedPaths != null) {
            preSelectPaths();
        }

    }

    //-----------------------------------------------------------------------------------------------------
    private void preSelectPaths() {
        String[] paths = preSelectedPaths.split(",");

    }

    //-----------------------------------------------------------------------------------------------------
    protected void initializeSelectionWidgetsAndVariables() {
        currentTreeSelection = new ArrayList<TreePath>();
        currentlySelectedExperiments = new ArrayList<MetaData>();
        loadSelectedConditionsButton.setText(" 0 ");
        loadSelectedConditionsButton.setEnabled(false);
        loadSelectedConditionsButton.setToolTipText("No conditions are currently selected");
        currentTreeSelection = new ArrayList<TreePath>();
        MetaData[] experiments = currentlySelectedExperiments.toArray(new MetaData[0]);
        for (MetaData experiment : experiments) experiment.clearSelectionCriteria();

    } // intializeSelectionWidgetsAndVariables

    //-----------------------------------------------------------------------------------------------------
    public void setMovieController(MovieController movieController) {
        this.movieController = movieController;
    }

    //-----------------------------------------------------------------------------------------------------
    public JFrame getMainFrame() {
        return jframe;
    }

    //-----------------------------------------------------------------------------------------------------
    public JToolBar getExtraToolBar() {
        return extraToolbar;
    }

    //-----------------------------------------------------------------------------------------------------
    public JToolBar getToolBar() {
        return toolbar;
    }

    //-----------------------------------------------------------------------------------------------------
    public JPanel getToolbarEtcPanel() {
        return toolbarEtcPanel;
    }

    //-----------------------------------------------------------------------------------------------------
    public JTree getJTree() {
        return treeWidget;
    }

    //-----------------------------------------------------------------------------------------------------
    public DataMatrixViewer getDataMatrixViewer() {
        return dataMatrixViewer;
    }

    //-----------------------------------------------------------------------------------------------------
    public int getSelectedConditionCount() {
        return selectedConditionCount;
    }

    //-----------------------------------------------------------------------------------------------------
    protected String createFrameTitle() {
        String id = "$Revision: 3545 $";
        String signature = "Revision: ";
        int start = id.indexOf(signature);
        start += signature.length();
        int end = id.indexOf(" $", start);
        String versionNumber = id.substring(start, end);
        return "Data Matrix Browser " + versionNumber;
    }

    //-----------------------------------------------------------------------------------------------------
    protected void createFrame() {
        jframe = new JFrame();
        jframe.setTitle(createFrameTitle());
        jframe.setJMenuBar(createMenuBar());
        jframe.setSize(800, 600);
        jframe.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        jframe.addWindowListener(this);
        MiscUtil.setApplicationIcon(jframe);


        jframe.getContentPane().add(createGui());

        // jframe.pack ();
        jframe.setVisible(true);
        ToolTipManager.sharedInstance().setInitialDelay(0);

    } // createFrame

    //-----------------------------------------------------------------------------------------------------
    protected JMenuBar createMenuBar() {
        menubar = new JMenuBar();
        JMenu fileMenu = new JMenu("File");
        fileMenu.setMnemonic('F');
        menubar.add(fileMenu);

        JMenuItem m0 = new JMenuItem("Load Tab-delimited matrix file...");
        m0.setMnemonic('L');
        m0.setAccelerator(
          KeyStroke.getKeyStroke(KeyEvent.VK_L, ActionEvent.META_MASK));
        JMenuItem m1 = new JMenuItem("Open New Repository (local file)");
        m1.setMnemonic('O');
        m1.setAccelerator(
          KeyStroke.getKeyStroke(KeyEvent.VK_O, ActionEvent.META_MASK));
        //JMenuItem m2 = new JMenuItem ("Save Current State...");
        JMenuItem m3 = new JMenuItem("Quit");
        m3.setMnemonic('Q');
        m3.setAccelerator(
          KeyStroke.getKeyStroke(KeyEvent.VK_Q, ActionEvent.META_MASK));



        fileMenu.add(m0);
        fileMenu.add(m1);
        fileMenu.add(new JMenuItem(new OpenRepositoryFromUrlAction()));
        fileMenu.add(new JMenuItem(new SaveCurrentlyOpenMatricesAction()));
        fileMenu.add(new JMenuItem(new LoadSavedMatrixProjectsAction()));
        fileMenu.add(m3);

        m0.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                chooseAndLoadNewDataFile();
            }
        });

        m1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                chooseNewRepositoryDirectory();
            }
        });

        //m2.addActionListener (new ActionListener () {
        //  public void actionPerformed (ActionEvent e) {
        //    saveAll ();
        //    }});

        m3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                exit();
            }
        });

        JMenu searchMenu = new JMenu("Search");
        searchMenu.setMnemonic('S');
        JMenuItem findItem = new JMenuItem("Find...");
        findItem.setMnemonic('F');
        findItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_F, ActionEvent.META_MASK));
        searchMenu.add(findItem);

        menubar.add(searchMenu);

        findItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                SearchDialog sd = new SearchDialog(repository);
                sd.pack();
                sd.setLocation(jframe.getLocation());
                sd.setVisible(true);
                if (sd.okPressed()) {
                    try {
                        MetaDataList mdl = new MetaDataList(repository);
                        Condition[] results = mdl.filterMetaDataByTags(sd.getTags());
                        if (results.length == 0) {
                            JOptionPane.showMessageDialog(jframe, "No results matched your query.");
                            actionPerformed(actionEvent);
                            return;
                        }
                        SearchResultsFrame resultsFrame = new SearchResultsFrame(sd.getTags(), results, TreeDataViewer.this);
                        resultsFrame.pack();
                        resultsFrame.setVisible(true);
                    } catch (Exception e) {
                        e.printStackTrace();
                        JOptionPane.showMessageDialog(jframe, "There was an error processing your search.");
                    }
                }
            }
        });


        return menubar;

    } // createMenuBar
//-----------------------------------------------------------------------------------------------------

    /**
     * expose ability to add menus, so we can add the Gaggle menu externally.
     * @param menu the menu to add
     */
    public void addMenu(JMenu menu) {
        menubar.add(menu);
    }

    //-----------------------------------------------------------------------------------------------------
    public void windowClosing(WindowEvent e) {
        exit();
    }

    //-----------------------------------------------------------------------------------------------------
    public void windowClosed(WindowEvent e) {
    }

    public void windowOpened(WindowEvent e) {
    }

    public void windowIconified(WindowEvent e) {
    }

    public void windowDeiconified(WindowEvent e) {
    }

    public void windowActivated(WindowEvent e) {
    }

    public void windowDeactivated(WindowEvent e) {
    }
//-----------------------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------------------
    protected void createTreeNodes(DefaultMutableTreeNode root, HashMap hash) {
        String[] kids = (String[]) hash.keySet().toArray(new String[0]);
        Arrays.sort(kids);

        for (String kid : kids) {
            DefaultMutableTreeNode newNode = new DefaultMutableTreeNode(kid);
            root.add(newNode);
            createTreeNodes(newNode, (HashMap) hash.get(kid));
        }

    } // createTreeNodes

    //-----------------------------------------------------------------------------------------------------
    protected void updateCurrentTreeSelection(TreeSelectionEvent treeSelectionEvent)
// update class member variable "currentTreeSelection" which holds a list
// of TreePaths (a TreePath is an array of Objects (here, Strings))
// specifying the tree branch.  paths are either added or removed
    {
        TreePath[] paths = treeSelectionEvent.getPaths();
        for (TreePath path : paths) {
            boolean added = treeSelectionEvent.isAddedPath(path);
            if (added)
                currentTreeSelection.add(path);
            else if (currentTreeSelection.contains(path))
                currentTreeSelection.remove(path);
        } // for p

    } // updateCurrentTreeSelection

    //-----------------------------------------------------------------------------------------------------
    protected HashMap<String, ArrayList<String>> convertTreePathsToHash(ArrayList<TreePath> currentTreeSelection)
// translate from, eg,
//
// [[Experiments, environmental, metals, FeSO4, time, -1],
//  [Experiments, environmental, metals, FeSO4, time, 20]]
//
// to
//
// {"environmental:metals:FeSO4": ["time:-1", "time:20"]}
//
    {
        HashMap<String, ArrayList<String>> result = new HashMap<String, ArrayList<String>>();
        if (experimentNavigator == null)
            return result;

        TreePath[] paths = currentTreeSelection.toArray(new TreePath[0]);
        //System.out.println ("convertTreePathsToHash, path count: " + paths.length);
        for (TreePath path : paths) {
            //System.out.println ("    path " + p + ": " + path);
            Object[] pathComponents = path.getPath();  // DefaultMutableTreeNode's
            // System.out.println ("pathComponents.length: " + pathComponents.length);
            String[] pathNames = new String[pathComponents.length - 1];
            //System.out.println ("   component count in path: " + pathComponents.length);
            for (int c = 1; c < pathComponents.length; c++) {
                //System.out.println ("    " + c + ": " + pathComponents [c] + " (" +
                //                    pathComponents[c].getClass ());
                pathNames[c - 1] = pathComponents[c].toString();
            }
            ArrayList experimentKeysList =
                    experimentNavigator.findExperimentKeyForPerturbation(pathNames);
            //System.out.println ("experiment keys, count: " + experimentKeysList.size ());

            for (Object anExperimentKeysList : experimentKeysList) {
                String[] tmp = (String[]) anExperimentKeysList;
                String experimentKey = tmp[0];
                String condition = tmp[1];
                if (condition == null)
                    condition = "";
                //System.out.println ("TDB.convertTreePathsToHash, experimentKey: " + experimentKey);
                //System.out.println ("TDB.convertTreePathsToHash,     condition: " + condition);
                if (!result.containsKey(experimentKey))
                    result.put(experimentKey, new ArrayList<String>());
                ArrayList<String> list = result.get(experimentKey);
                if (!list.contains(condition))
                    list.add(condition);
                result.put(experimentKey, list);
            } // for k
        } // for p

        return result;

    } // convertTreePathsToHash

    //-----------------------------------------------------------------------------------------------------
    public void valueChanged(TreeSelectionEvent treeSelectionEvent)
// respond to a change in the selection state of the treeWidget.
// the goal is to translate the current selection state of the JTree widget
// into selected conditions in all of the relevant experiements -- which are
// here actually MetaData objects.  these MetaData objects originally supplied
// the JTree hierarchy.  they also carry information about every observed condition
// in the experiment (typically a dosage level, or a point in a time course);
// any of these conditions may be selected, and that's the result we seek here:
// that all JTree selections are mirrored by selections in the MetaData objects.
// these MetaData objects can then be used to create matrices that consist only
// of the selected columns (data columns are, in effect, experimental conditions)
    {
        updateCurrentTreeSelection(treeSelectionEvent);

        if (currentTreeSelection.size() == 0) {
            initializeSelectionWidgetsAndVariables();
            MetaData[] experiments = currentlySelectedExperiments.toArray(new MetaData[0]);
            for (MetaData experiment : experiments) experiment.clearSelectionCriteria();
            return;
        }

        // there is at least one selection.  translate the TreePaths into a hash.

        loadSelectedConditionsButton.setEnabled(true);
        experimentDetailsAction.setEnabled(true);
        HashMap<String, ArrayList<String>> selectedConditions = convertTreePathsToHash(currentTreeSelection);

        if (selectedConditions == null || selectedConditions.size() == 0) {
            loadSelectedConditionsButton.setText(" 0 ");
            loadSelectedConditionsButton.setEnabled(false);
            loadSelectedConditionsButton.setToolTipText("No conditions are currently selected");

            experimentDetailsAction.setEnabled(false);

            currentTreeSelection = new ArrayList<TreePath>();
            MetaData[] experiments = currentlySelectedExperiments.toArray(new MetaData[0]);
            for (MetaData experiment : experiments) experiment.clearSelectionCriteria();
            currentlySelectedExperiments.clear();
            return;
        }

        currentlySelectedExperiments = new ArrayList<MetaData>();
        selectedConditionCount = 0;

        String[] experimentKeys = selectedConditions.keySet().toArray(new String[0]);
        MetaData[] experimentList = new MetaData[0];
        Map allAliases = new HashMap();

        for (String experimentKey : experimentKeys) {
            //System.out.println ("  e = " + e + "  " + experimentKey);
            // the experiment key is just a convenenient form of the tree paths selected in
            // the jtree. it might be an explicit, single experiment
            //    "environmental:metals:FeSO4"
            // or implicitly, a whole group of experiments
            //    "environmental:metals"
            // which implies FeSO4, cobalt, copper, iron, manganese, ...

            experimentList = experimentNavigator.getExperimentByPerturbationList(experimentKey);

            // now select conditions in each experiment
            for (MetaData experiment : experimentList) {
                experiment.clearSelectionCriteria();
                if (!currentlySelectedExperiments.contains(experiment))
                    currentlySelectedExperiments.add(experiment);
                ArrayList<String> conditions = selectedConditions.get(experimentKey);
                experiment = makeSelections(experiment, conditions);
                for (String alias : experiment.getSelectedConditionsAsAliases()) {
                    allAliases.put(alias, "-");
                }
            } // for i
        } // for e
        selectedConditionCount = allAliases.keySet().size();

        if (experimentKeys.length == 1) { // we can only read more about one experiment at a time
            followExperimentLinksButton.setEnabled(true);
            MetaData experimentMetaData = experimentList[0];
            selectedExperimentInfoUrls = experimentMetaData.getReferenceLinks();
        } else
            followExperimentLinksButton.setEnabled(false);



        String tmp = String.valueOf(selectedConditionCount);
        loadSelectedConditionsButton.setText(tmp);
        loadSelectedConditionsButton.setToolTipText("Load the " + tmp + " currently selected conditions");

    } // valueChanged


    //-----------------------------------------------------------------------------------------------------
    protected MetaData makeSelections(MetaData experiment, ArrayList conditions) {

// selected conditions for the current experiment take several forms,
// directly reflecting how deeply into the JTree the user clicked
//   treePath                             experimentKey                 conditions
//  --------------------------          -------------------------       ----------
// environmental:metals:FeSO4           environmental:metals:FeSO4         null
// environmental:metals:FeSO4:time      environmental:metals:FeSO4         time
// environmental:metals:FeSO4:time:-1   environmental:metals:FeSO4         time:-1

        //System.out.println ("TDB.makeSelections, conditions count: " + conditions.size ());
        //System.out.println ("type of conditions.get (0): " + conditions.get(0).getClass ());
        String[] conditionRawPair = (String[]) conditions.toArray(new String[0]);
        //for (int i=0; i < conditionRawPair.length; i++)
        //System.out.println ("conditionRawPair " + i + ": " + conditionRawPair [i]);

        if (conditionRawPair.length == 0)
            experiment.selectAllConditions();
        else if (conditionRawPair.length == 1 && conditionRawPair[0].length() == 0)
            experiment.selectAllConditions();
        else {
            for (String aConditionRawPair : conditionRawPair) {
                String[] conditionPair = aConditionRawPair.split(":");
                if (conditionPair.length == 1)   // select all values for this condition
                    experiment.selectConditionByName(conditionPair[0]);
                else if (conditionPair.length == 2) {
                    String name = conditionPair[0];
                    String value = conditionPair[1];
                    experiment.addSelectionCriterion(name, value);
                }
            } // for i
        }
        //System.out.println ("TDV.makeSelections, species: " + experiment.getSpecies ());
        return experiment;

    } // makeSelections

    //-----------------------------------------------------------------------------------------------------
    protected void expandOrCollapseTree(boolean expand) {
        boolean done = false;
        while (!done) {
            int rowCount = treeWidget.getRowCount();
            for (int r = rowCount; r > 0; r--) {
                if (expand) {
                    treeWidget.expandRow(r);
                } else {
                    treeWidget.collapseRow(r);
                }
            }
            if (rowCount == treeWidget.getRowCount())
                done = true;
        } // while !done

    } // expandTree

    //-----------------------------------------------------------------------------------------------------
    public void addRepository(String repository) {
        DefaultMutableTreeNode root = new DefaultMutableTreeNode("Experiments");

        if (repository != null && repository.length() > 0) try {
            System.out.println("repository: " + repository);
            experimentNavigator = new MetaDataNavigator(repository);
            HashMap experimentsTree = experimentNavigator.getTree();
            createTreeNodes(root, experimentsTree);
        }
        catch (Exception ex0) {
            ex0.printStackTrace();
            JOptionPane.showMessageDialog(jframe, ex0.getMessage(), "Add Repository Error",
                    JOptionPane.ERROR_MESSAGE);
        } // catch

        treeWidget = new JTree(root);
        treeWidget.getSelectionModel().setSelectionMode(TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION);
        treeWidget.addTreeSelectionListener(this);

        treePane.setViewportView(treeWidget);
        initializeSelectionWidgetsAndVariables();

    } // addRepository

    //-----------------------------------------------------------------------------------------------------
    protected JPanel createGui() {
        JPanel panel = new JPanel();
        panel.setMinimumSize(new Dimension(600, 400));

        panel.setLayout(new BorderLayout());

        JPanel chooserPanel = new JPanel();
        chooserPanel.setLayout(new BorderLayout());
        treePane = new JScrollPane();

        JPanel treePanel = new JPanel();
        treePanel.setLayout(new BorderLayout());
        treePanel.setMinimumSize(new Dimension(200, 10));
        treePanel.add(treePane, BorderLayout.CENTER);

        JPanel selectionPanel = new JPanel();
        //JLabel label = new JLabel ("Selection: ");
        //selectionPanel.add (label);
        conditionCounterTextField = new JTextField(5);
        conditionCounterTextField.setToolTipText("Number of columns explicitly selected in tree");
        selectionPanel.add(conditionCounterTextField);
        //treePanel.add (selectionPanel, BorderLayout.SOUTH);

        resultsPanel = new JPanel();
        resultsPanel.setMinimumSize(new Dimension(600, 400));
        resultsPanel.setLayout(new BorderLayout());

        outerSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, treePanel, resultsPanel);
        outerSplitPane.setDividerLocation(200);
        outerSplitPane.setOneTouchExpandable(true);
        dataBrowserHasEntireWindow = false;

        panel.add(outerSplitPane, BorderLayout.CENTER);

        toolbar = new JToolBar();
        extraToolbar = new JToolBar();
        JPanel placementPanel1 = new JPanel();
        JPanel placementPanel2 = new JPanel();
        placementPanel1.setLayout(new BorderLayout());
        placementPanel2.setLayout(new BorderLayout());
        placementPanel2.add(toolbar, BorderLayout.WEST);
        placementPanel1.add(extraToolbar, BorderLayout.WEST);
        //placementPanel1.setBackground (Color.BLUE);
        //placementPanel2.setBackground (Color.GREEN);
        JPanel toolbarsPanel = new JPanel();
        toolbarsPanel.setLayout(new GridLayout(2, 1));
        toolbarsPanel.add(placementPanel1);
        toolbarsPanel.add(placementPanel2);

        toolbar.setFloatable(false);
        toolbarEtcPanel = new JPanel();
        toolbarEtcPanel.setLayout(new BorderLayout());
        toolbarEtcPanel.add(toolbarsPanel, BorderLayout.WEST);
        panel.add(toolbarEtcPanel, BorderLayout.NORTH);
        //toolbar.add (JButton ("Add Data Source", actionPerformed=addSource))

        newRepositoryButton = new JButton("New...");
        newRepositoryButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                chooseNewRepositoryDirectory();
            }
        });
        //toolbar.add (newRepositoryButton);

        saveAllButton = new JButton("Save...");
        saveAllButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                saveAll();
            }
        });
        //toolbar.add (saveAllButton);

        loadSelectedConditionsButton = new JButton(" 0 "); // IconFactory.getLoadSelectedConditionsIcon());
        loadSelectedConditionsButton.setBackground(Color.WHITE);
        loadSelectedConditionsButton.setForeground(Color.RED);
        loadSelectedConditionsButton.setToolTipText("No conditions are currently selected");
        loadSelectedConditionsButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                loadSelectedConditions();
            }
        });

        loadSelectedConditionsButton.setBackground(Color.WHITE);
        loadSelectedConditionsButton.setEnabled(false);

        JButton expandButton = new JButton("+/-"); // IconFactory.getExpandAndContractIcon());
        expandButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                expandContractDataBrowser();
            }
        });


        expandButton.setToolTipText("Toggle Tree Visibility");
        expandButton.setBackground(Color.WHITE);
        toolbar.add(expandButton);

        JButton expandTreeButton = new JButton("T");
        expandTreeButton.setToolTipText("Expand/Collapse Tree");
        expandTreeButton.setBackground(Color.WHITE);
        expandTreeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent actionEvent) {
                expandOrCollapseTree(treeCollapsed);
                treeCollapsed = !treeCollapsed;
            }
        } );
        toolbar.add(expandTreeButton);

        showHideWorkFramesButton = new JButton("H/S");
        showHideWorkFramesButton.setToolTipText("Hide/Show Data");

        showHideWorkFramesButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                toggleVisibilityOfWorkFrames();
            }
        });
        toolbar.add(showHideWorkFramesButton);

        JButton dismissButton = new JButton("Quit");
        dismissButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                exit();
            }
        });
        dismissButton.setBackground(Color.WHITE);
        dismissButton.setToolTipText("Close This Window");
        // toolbar.add (dismissButton);

        toolbar.add(loadSelectedConditionsButton);

        followExperimentLinksButton = new JButton("WWW");
        followExperimentLinksButton.setEnabled(false);
        String msg = "<html>Experiment information -- <br> select just one experiment</html>";
        followExperimentLinksButton.setToolTipText(msg);
        followExperimentLinksButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                followExperimentWebLinks();
            }
        });

        toolbar.add(followExperimentLinksButton);

        experimentDetailsAction = new ExperimentDetailsAction();
        experimentDetailsAction.setEnabled(false);
        JButton experimentDetailsButton = new JButton(experimentDetailsAction);
        toolbar.add(experimentDetailsButton);

        //todo deal with enabling/disabling
        closeAllTabsAction = new CloseAllTabsAction();
        //closeAllTabsAction.setEnabled(false);
        JButton closeAllTabsButton = new JButton(closeAllTabsAction);
        toolbar.add(closeAllTabsButton);

        // TODO move this? Not really part of creating gui
        addRepository(repository);

        return panel;

    } // createGui
//-----------------------------------------------------------------------------------------------------
    //todo? change so that you can see a list of existing projects (to overwrite one), OR enter a new name
    class SaveCurrentlyOpenMatricesAction extends AbstractAction {
        SaveCurrentlyOpenMatricesAction() {
            super("Save currently open matrices...");
            putValue(ACCELERATOR_KEY, KeyStroke.getKeyStroke(KeyEvent.VK_S, ActionEvent.META_MASK));
            int ch = 'S';
            putValue(MNEMONIC_KEY, ch);
            this.putValue(Action.SHORT_DESCRIPTION, "Save all currently open matrices to reopen later.");
        }


        public void actionPerformed(ActionEvent e) {
            if (dataMatrixViewer.getAllViews().length == 0) {
                JOptionPane.showMessageDialog(jframe, "No tabs to save!");
                return;
            }

            SaveMatrixProjectDialog dialog = new SaveMatrixProjectDialog();
            dialog.pack();
            dialog.setVisible(true);

            if (!dialog.isOk()) {
                return;
            }

            final MatrixPersister mp = new MatrixPersister();
            MatrixProject already = mp.getProjectByName(dialog.getName());
            if (already != null) { // todo test this
                JOptionPane.showMessageDialog(jframe, "Name is already taken; choose another.");
                dialog.setVisible(true);
                actionPerformed(null);
                return;
            }
            mp.closeDatabase();


            final MatrixProject project = new MatrixProject();
            project.setName(dialog.getName());
            project.setDescription(dialog.getDescription());
            project.setDate(new Date());
            project.setViews(dataMatrixViewer.getAllMatrixSpreadsheets());

            SwingUtilities.invokeLater(new Runnable(){
                public void run() {
                    try {
                        mp.saveProject(project);
                    } catch (Exception e) {
                        e.printStackTrace();
                        JOptionPane.showMessageDialog(jframe,"There was a problem saving your data.\n"
                        + "It may or may not be saved.\n"
                        + "Go to the 'Load Saved Matrix Project' menu item\n"
                        + "to check.");
                    }
                }
            });
        }
    }

    class LoadSavedMatrixProjectsAction extends AbstractAction {
        LoadSavedMatrixProjectsAction() {
            super("Load saved matrix project...");
            putValue(ACCELERATOR_KEY, KeyStroke.getKeyStroke(KeyEvent.VK_P, ActionEvent.META_MASK));
            int ch = 'P';
            putValue(MNEMONIC_KEY, ch);
            this.putValue(Action.SHORT_DESCRIPTION, "Load a previously saved set of matrices.");
            
        }

        public void actionPerformed(ActionEvent e) {
            MatrixPersister mp = new MatrixPersister();
            if (!mp.doesDatabaseExist()) {
                JOptionPane.showMessageDialog(jframe, "No saved projects to load!");
                return;
            }
            MatrixProject[] allProjects = mp.getAllProjects();
            //MatrixProject[] allProjects = perst.toArray(new MatrixProject[0]);


            if (allProjects == null || allProjects.length == 0) {
                JOptionPane.showMessageDialog(jframe, "No saved projects to load!");
                return;
            }
            LoadMatrixProjectDialog dialog = new LoadMatrixProjectDialog(allProjects);
            dialog.pack();
            dialog.setVisible(true);

            if (dialog.isOk()) {
                if (dialog.closeAllTabsFirst()) {
                    closeAllTabsAction.actionPerformed(null);
                    // todo - determine whether we want to do the following whether the checkbox was checked or not
                    currentlySelectedExperiments = new ArrayList<MetaData>();
                    currentTreeSelection = new ArrayList<TreePath>();
                    treeWidget.clearSelection();
                    initializeSelectionWidgetsAndVariables();
                    expandOrCollapseTree(false);
                }
                for (DataMatrixView view : dialog.getSelectedProject().getViews()) {
                    dataMatrixViewer.addView(view);
                }
                jframe.pack();
                outerSplitPane.setDividerLocation(200);
                // todo - clear selected conditions from tree
            }
            mp.closeDatabase();
        }
    }

    class CloseAllTabsAction extends AbstractAction {
        public CloseAllTabsAction() {
            super("X");
            this.putValue(Action.SHORT_DESCRIPTION, "Close all open tabs");
        }

        public void actionPerformed(ActionEvent e) {
            if (dataMatrixViewer.getTabCount() == 0) {
                JOptionPane.showMessageDialog(jframe, "No Tabs to close!");
                return;
            }
            for (int i = (dataMatrixViewer.getTabCount()-1); i > -1; i--) {
                dataMatrixViewer.closeTab(i);
                System.gc ();
                dataMatrixViewer.updateCongruentMatrixInfo(i-1);
            }
        }
    }

    /**
     * show details (metadata) about selected experiment(s).
     */
    class ExperimentDetailsAction extends AbstractAction {

        ExperimentDetailsAction() {
            super("Details");
            this.putValue(Action.SHORT_DESCRIPTION, "View details for selected experiment(s).");
        }

        public void actionPerformed(ActionEvent e) {
            MetaDataViewer mdv = new MetaDataViewer(jframe);
            mdv.setExperimentList(currentlySelectedExperiments);
            mdv.setVisible(true);
        }
    }

    //-----------------------------------------------------------------------------------------------------
    public void followExperimentWebLinks() {
        String[] linkTypes = (String[]) selectedExperimentInfoUrls.keySet().toArray(new String[0]);

        for (String linkType : linkTypes) {
            String url = (String) selectedExperimentInfoUrls.get(linkType);
            MiscUtil.displayWebPage(url);
        }


    } // followExperimentWebLinks

    //-----------------------------------------------------------------------------------------------------
    public void toggleVisibilityOfWorkFrames() {
        workFramesVisible = !workFramesVisible;
        //String label = "-";
        //if (!workFramesVisible)
        // label = "+";
        //showHideWorkFramesButton.setText (label);
        outerSplitPane.setVisible(workFramesVisible);
        outerSplitPane.setDividerLocation(200);
        jframe.pack();
    }

    //---------------------------------------------------------------------------------
    protected void chooseAndLoadNewDataFile() {
        JFileChooser chooser = new JFileChooser(settings.getWorkingDirectory());
        // chooser.setFileFilter (new RepositoryFilter ());
        chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
        int returnVal = chooser.showOpenDialog(jframe);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File chosenFile = chooser.getSelectedFile();
            if (chosenFile != null) {
                settings.setWorkingDirectory(chooser.getCurrentDirectory().getPath());
                String speciesForNewData = getSpeciesName(getMainFrame());
                try {
                    addMatrix(chosenFile, speciesForNewData);
                }
                catch (Exception ex0) {
                    ex0.printStackTrace();
                    JOptionPane.showMessageDialog(jframe, "Read Matrix Error: " + ex0.getMessage(), ex0.getMessage(),
                            JOptionPane.ERROR_MESSAGE);
                }
            }
        }
    }

    //---------------------------------------------------------------------------------
/**
 * put up possibly two dialog boxes.  the first presents known organisms (unfortunately
 * hard-coded; must improve upon that!), 'unknown' and 'other...'.  choosing the last
 * leads to the second dialog box, which prompts for the organims name using free text
 * entry -- no sanity checking provided for this option.
 * <p/>
 * @param frame the frame
 * @return species name
 * todo: obtain known species from some configured source.  in the absence of a
 * todo: configured source, present just "unknown" and "other..."
 */
    protected String getSpeciesName(JFrame frame) {
        String result = "unknown";

        String[] knownSpecies = {"Halobacterium sp.",
                "Helicobacter pylori",
                "Homo sapiens",
                "Mus musculus",
                "Drosophila melanogaster",
                "unknown",
                "other..."};

        String choice = (String) JOptionPane.showInputDialog(
                frame, "Specify the organism for these data:",
                "Choose Species",
                JOptionPane.PLAIN_MESSAGE,
                null,
                knownSpecies,
                "unknown");

        if (choice == null || choice.equals("unknown"))
            result = "unknown";
        else if (!choice.equals("other..."))
            result = choice;
        else {
            String newSpecies = (String) JOptionPane.showInputDialog(
                    frame, "Specify the organism for these data:",
                    "Enter species name",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    null,
                    "");
            if (newSpecies != null && newSpecies.length() > 0)
                result = newSpecies;
        } // if "other..."

        return result;

    }  // getSpeciesName

    //-----------------------------------------------------------------------------------
    protected void addMatrix(File file, String speciesForNewData) throws Exception {
        DataMatrixFileReader reader = new DataMatrixFileReader(file.getPath());
        reader.read();
        org.systemsbiology.gaggle.core.datatypes.DataMatrix matrix = reader.get();
        matrix.setSpecies(speciesForNewData);
        dataMatrixViewer.addMatrixSpreadsheetView(matrix, experimentNavigator);

    } // addMatrix

    //---------------------------------------------------------------------------------
    protected void chooseNewRepositoryDirectory() {
        currentTreeSelection = new ArrayList<TreePath>();

        JFileChooser chooser = new JFileChooser(settings.getWorkingDirectory());
        chooser.setFileFilter(new RepositoryFilter());
        chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        int returnVal = chooser.showOpenDialog(jframe);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File chosenDirectory = chooser.getSelectedFile();
            if (chosenDirectory != null) {
                // System.out.println ("You chose directory: " + chosenDirectory);
                settings.setWorkingDirectory(chooser.getCurrentDirectory().getPath());
                //System.out.println ("currentDirectory: " + currentDirectory);
                addRepository(chosenDirectory.getPath());
            } // if !null
        } // if approve

    } // chooseNewRepositoryDirectory

    //-----------------------------------------------------------------------------------------------------
    protected void saveAll() {
        JFileChooser chooser = new JFileChooser(settings.getWorkingDirectory());
        // chooser.setFileFilter (new RepositoryFilter ());
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        int returnVal = chooser.showSaveDialog(jframe);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File chosenDirectory = chooser.getSelectedFile();
            if (chosenDirectory != null) {
                settings.setWorkingDirectory(chooser.getCurrentDirectory().getPath());
                //System.out.println ("save to directory: " + chosenDirectory.getAbsolutePath ());
                if (!chosenDirectory.exists())
                    try {
                        System.out.println("about to mkdir: " + chosenDirectory);
                        chosenDirectory.mkdir();
                        System.out.println("after mkdir");
                    }
                    catch (Exception ex0) {
                        ex0.printStackTrace();
                        JOptionPane.showMessageDialog(jframe, "Save Error", ex0.getMessage(),
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                System.out.println("  exists??? " + chosenDirectory.exists());
                dataMatrixViewer.saveAll(chosenDirectory);
            } // if choseDirectory != null
        } // if approve

    } // saveAll

    //-----------------------------------------------------------------------------------------------------
    class RepositoryFilter extends FileFilter {
        public boolean accept(File f) {
            if (f.isDirectory()) return true;
            String name = f.getName();
            return name.endsWith(".xml");
        }

        public String getDescription() {
            return "directories and xml files";
        }
    }

    //-----------------------------------------------------------------------------------------------------
    public void exit() {
        notifyExitListeners();
        dataMatrixViewer = null;
        jframe.dispose();
        System.exit(0);

    }

    //---------------------------------------------------------------------------------
/**
 * listeners will be notified when the user stops the program
 * @param listener the listener
 */
    public void addExitListener(ExitListener listener) {
        synchronized (exitListeners) { // synchronization here unlikely to have useful semantics
            exitListeners.add(listener);
        }
    }

    //---------------------------------------------------------------------------------
    public void removeExitListener(ExitListener listener) {
        synchronized (exitListeners) { // synchronization here unlikely to have useful semantics
            exitListeners.remove(listener);
        }
    }

    //---------------------------------------------------------------------------------
    public void notifyExitListeners() {
        synchronized (exitListeners) { // synchronization here unlikely to have useful semantics
            for (ExitListener listener : exitListeners) {
                listener.exitInProgress();
            }
        }
    }
//-----------------------------------------------------------------------------------------------------
/*
public class AddRepositoryAction extends AbstractAction {

  AddRepositoryAction () {
    super ("Add Repository");
  }

  public void actionPerformed (ActionEvent e) {
    JFileChooser chooser = new JFileChooser (settings.getWorkingDirectory());
    chooser.setFileSelectionMode (JFileChooser.DIRECTORIES_ONLY);
    int status = chooser.showOpenDialog (jframe);
    if (status == JFileChooser.APPROVE_OPTION) {
      settings.setWorkingDirectory( chooser.getCurrentDirectory().getPath() );
      // wasTODO
    }
  }

}
*/
//-----------------------------------------------------------------------------------------------------

    /**
     * Allow user to load a remote repository using either http:// or httpIndirect://
     * protocols.
     */
    class OpenRepositoryFromUrlAction extends AbstractAction {

        public OpenRepositoryFromUrlAction() {
            super("Open New Repository (URL)");
            putValue(ACCELERATOR_KEY, KeyStroke.getKeyStroke(KeyEvent.VK_U, ActionEvent.META_MASK));
            int ch = 'U';
            putValue(MNEMONIC_KEY, ch);
            this.putValue(Action.SHORT_DESCRIPTION, "Open a new repository located by URL.");
        }

        public void actionPerformed(ActionEvent e) {
            String url = JOptionPane.showInputDialog(jframe,
                    "Repository URL:", "Enter Repository URL", JOptionPane.QUESTION_MESSAGE);
            if (url != null) {
                addRepository(url);
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------
    public void loadSelectedConditions() {
        org.systemsbiology.gaggle.core.datatypes.DataMatrix[] matrices = new org.systemsbiology.gaggle.core.datatypes.DataMatrix[0];
        try {
            matrices = combineSelectedConditions();
        }
        catch (Exception ex0) {
            ex0.printStackTrace();
            String msg = "<html>Error combining selected conditions: <br>" + ex0.getMessage() + "</html>";
            JOptionPane.showMessageDialog(jframe, msg, "Data Error", JOptionPane.ERROR_MESSAGE);
        }

        putMatricesIntoDataMatrixViewer(matrices);

        dataMatrixViewer.addMetaData(currentlySelectedExperiments);
        if (movieController != null) {
            movieController.loadMatrices();
        }
        jframe.pack();
        outerSplitPane.setDividerLocation(200);

    } // loadSelectedConditions

    private void putMatricesIntoDataMatrixViewer(org.systemsbiology.gaggle.core.datatypes.DataMatrix[] matrices) {
        if (matrices.length == 0)
            return;
        if (dataMatrixViewer != null) {
            for (org.systemsbiology.gaggle.core.datatypes.DataMatrix matrix : matrices) dataMatrixViewer.addMatrixSpreadsheetView(matrix, experimentNavigator);
        }
        else {
            try {
                dataMatrixViewer = new DataMatrixViewer(matrices, experimentNavigator);
                resultsPanel.add(dataMatrixViewer, BorderLayout.CENTER);
                resultsPanel.paintAll(resultsPanel.getGraphics());
            }
            catch (Exception ex0) {
                ex0.printStackTrace();
            }
        } // else:  create a new DataMatrixViewer
        return;
    }

    //-----------------------------------------------------------------------------------------------------
    public void expandContractDataBrowser() {
        if (dataBrowserHasEntireWindow) {
            int lastLocation = outerSplitPane.getLastDividerLocation();
            outerSplitPane.setDividerLocation(lastLocation);
            dataBrowserHasEntireWindow = false;
        } else {
            outerSplitPane.setDividerLocation(0);
            dataBrowserHasEntireWindow = true;
        }

    } // expandContractDataBrowser

    //-----------------------------------------------------------------------------------------------------
    protected org.systemsbiology.gaggle.core.datatypes.DataMatrix[] combineSelectedConditions() throws Exception
// create a single matrix from all of the selected columns in all selected experiments
    {
        if (currentlySelectedExperiments.size() == 0) {
            return new org.systemsbiology.gaggle.core.datatypes.DataMatrix[0];
        }

        HashMap<String, ArrayList> allSelectedMatrices = new HashMap<String, ArrayList>();


        for (MetaData experiment : currentlySelectedExperiments) {
            String[] columnNames = experiment.getSelectedConditionsAsAliases();
            if (columnNames.length > 0) {
                MatrixSlicer slicer = new MatrixSlicer(experiment);
                //try {
                HashMap slicedMatrices = slicer.slice();  // there may be one, two, or more...
                String[] matrixTypes = (String[]) slicedMatrices.keySet().toArray(new String[0]);
                for (String matrixType : matrixTypes) {
                    org.systemsbiology.gaggle.core.datatypes.DataMatrix m = (org.systemsbiology.gaggle.core.datatypes.DataMatrix) slicedMatrices.get(matrixType);
                    if (!allSelectedMatrices.containsKey(matrixType))
                        allSelectedMatrices.put(matrixType, new ArrayList());
                    (allSelectedMatrices.get(matrixType)).add(m);
                }
            }
        }

        ArrayList<org.systemsbiology.gaggle.core.datatypes.DataMatrix> finalMatrices = new ArrayList<org.systemsbiology.gaggle.core.datatypes.DataMatrix>();
        String[] matrixTypes = allSelectedMatrices.keySet().toArray(new String[0]);
        for (String matrixType : matrixTypes) {
            ArrayList matrixListOfOneType = allSelectedMatrices.get(matrixType);
            org.systemsbiology.gaggle.core.datatypes.DataMatrix[] matricesOfOneType = (org.systemsbiology.gaggle.core.datatypes.DataMatrix[]) matrixListOfOneType.toArray(new org.systemsbiology.gaggle.core.datatypes.DataMatrix[0]);
            MatrixCombiner combiner = new MatrixCombiner(matricesOfOneType);
            org.systemsbiology.gaggle.core.datatypes.DataMatrix combinedMatrix = combiner.combine();
            combinedMatrix.setFullName("(from TreeDataViewer)");
            combinedMatrix.setShortName(matrixType);
            combinedMatrix.setDataTypeBriefName(matrixType);
            finalMatrices.add(combinedMatrix);
        }

        return finalMatrices.toArray(new org.systemsbiology.gaggle.core.datatypes.DataMatrix[0]);

    } // combineSelectedConditions

    //---------------------------------------------------------------------------------
    public static void main(String[] args) {
        String repositoryName = null;
        if (args.length == 1)
            repositoryName = args[0];

        /*TreeDataViewer tdb =*/
        new TreeDataViewer(repositoryName);


    }
//---------------------------------------------------------------------------------
} // class TreeDataViewer
