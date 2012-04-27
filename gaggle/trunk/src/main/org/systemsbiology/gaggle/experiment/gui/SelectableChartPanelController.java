// SelectableChartPanelController.java
//-----------------------------------------------------------------------------------------------
/*
 * Copyright (C) 2006 by Institute for Systems Biology,
 * Seattle, Washington, USA.  All rights reserved.
 *
 * This source code is distributed under the GNU Lesser
 * General Public License, the text of which is available at:
 *   http://www.gnu.org/copyleft/lesser.html
 */

package org.systemsbiology.gaggle.experiment.gui;

import java.awt.Rectangle;

public interface SelectableChartPanelController {
    
    public void setDraggedRectangle(Rectangle draggedRectangle); 
}
