# makefile for csplugins/isb/pshannon/experiment/gui/
#-------------------------------------------------------------------------
# $Revision: 375 $
# $Date: 2005/04/13 16:53:01 $
# $Author: pshannon $
#-------------------------------------------------------------------------
JC=javac 
JI=java

.SUFFIXES: .java .class
.java.class: 
	$(JC) $<

default: objs subdirs

OBJS=   CheckBoxHeader.class \
        DataMatrixTableModel.class \
        MatrixViewCoordinator.class \
        DataMatrixView.class \
        MatrixViewCoordinator.class \
        MatrixSpreadsheet.class \
        DataMatrixViewer.class \
        TreeDataViewer.class \
        MutableList.class \
        SelectableChartPanel.class \
        SelectableChartPanelController.class \
        SpreadsheetLineFinder.class


objs: $(OBJS)

subdirs:
	(cd movie; make)
	(cd actions; make)
	(cd movie; make)

clean:
	- find . -name "*.class" -exec rm  {} ';'

REPOS=httpIndirect://db.systemsbiology.net:8080/halo/DataFetcher.py

# --- tree
# run the dmv with the experiemnt chooser tree on the left.  this is the gaggle's full
# DMV, lacking only the gaggle connection
#

tree: $(OBJS)
	java -Xmx1G org.systemsbiology.gaggle.experiment.gui.TreeDataViewer $(REPOS)




help:
	egrep "^#" makefile | sed "s/^#//"





