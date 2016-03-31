package plugins.danyfel80.segmentation.powerwatershed;

import icy.gui.dialog.MessageDialog;
import icy.roi.ROI;
import icy.roi.ROI2D;
import icy.sequence.Sequence;

import java.util.List;

import algorithms.danyfel80.segmentation.SegmentationAlgorithm;
import algorithms.danyfel80.segmentation.graphcut.GraphCutSegmentation;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.ylemontag.histogram.BadHistogramParameters;

/**
 * @author Daniel Felipe Gonzalez Obando
 * 
 *         Power Watershed Framework Plugin
 * 
 *         This plugin implements the segmentation framework proposed in the
 *         article
 *         "Power Watershed: A Unifying Graph-Based Optimization Framework" by
 *         Camille Couprie, Leo Grady, Laurent Najman, Hugues Talbot.
 */
public class PowerWatershedPlugin extends EzPlug {

	private EzVarSequence inSequence;
	private EzVarSequence inSeedsSequence;
	private EzVarInteger inP;
	private EzVarBoolean inIsPInfinite;
	private EzVarInteger inQ;
	private EzVarBoolean inIsQInfinite;
	private EzVarBoolean inUse8Connected;

	/* (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#initialize()
	 */
	@Override
	protected void initialize() {
		
		// Input instantiation and setup
		inSequence = new EzVarSequence("Sequence");
		inSequence.setOptional(false);
		inSequence.setToolTipText("The sequence to be segmented.");

		inSeedsSequence = new EzVarSequence("Seeds");
		inSeedsSequence.setOptional(true);
		inSeedsSequence.setToolTipText("The seeds for the segmentation, "
				+ "at least two classes must be specified.");

		inP = new EzVarInteger("P", 0, 0, 10, 1);
		inP.setOptional(false);
		inP.setToolTipText("The power of weights in the graph.");

		inIsPInfinite = new EzVarBoolean("P = inifity", false);
		inIsPInfinite.setOptional(false);
		inIsPInfinite.setToolTipText("If true the value of P is considered as "
				+ "infinite.");
		inIsPInfinite.addVarChangeListener(new EzVarListener<Boolean>() {
			@Override
			public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
				if (newValue) {
					inP.setValue(-1);
					inP.setVisible(false);
				} else {
					inP.setValue(0);
					inP.setVisible(true);
				}
			}
		});

		inQ = new EzVarInteger("Q", 1, 1, 2, 1);
		inQ.setOptional(false);
		inQ.setToolTipText("The power of neighboring differences in the graph.");

		inIsQInfinite = new EzVarBoolean("Q = inifity", false);
		inIsQInfinite.setOptional(false);
		inIsQInfinite.setToolTipText("If true the value of Q is considered as "
				+ "infinite.");
		inIsQInfinite.addVarChangeListener(new EzVarListener<Boolean>() {
			@Override
			public void variableChanged(EzVar<Boolean> source, Boolean newValue) {
				if (newValue) {
					inQ.setValue(0);
					inQ.setVisible(false);
				} else {
					inQ.setValue(1);
					inQ.setVisible(true);
				}
			}
		});
		
		inUse8Connected = new EzVarBoolean("Use 8 connected image", false);

		// Input components addition to UI
		addEzComponent(inSequence);
		addEzComponent(inSeedsSequence);
		addEzComponent(inP);
		addEzComponent(inIsPInfinite);
		addEzComponent(inQ);
		addEzComponent(inIsQInfinite);
		addEzComponent(inUse8Connected);
	}

	/* (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#execute()
	 */
	@Override
	protected void execute() {
		
		// Input validation
		if (validateInput() != 0) {
			return;
		}
		
		// Get the used algorithm
		SegmentationType algoType = SegmentationType.getAlgorithm(inP.getValue(), inQ.getValue());
		SegmentationAlgorithm algo = null;
		
		long startTime = 0, endTime = 0;
		
		switch (algoType) {
		case CollapseToSeeds:
			
			break;
		case l2NormVoronoi:
			
			break;
		case l1NormVoronoi:
	
			break;
		case GraphCuts:
			double lambda = 1;
			startTime = System.nanoTime();
			algo = new GraphCutSegmentation(inSequence.getValue(), lambda, inUse8Connected.getValue());
			endTime = System.nanoTime();
			break;
		case RandomWalker:
			
			break;
		case PowerWatershedQ1:
			
			break;
		case PowerWatershedQ2:
			
			break;
		case ShortestPathForest:
			
			break;
		}
		
		System.out.println("Time to prepare algorithm: " + (endTime - startTime)/1000000 + "ms");
		
		// Execute algorithm
		try {
		  startTime = System.nanoTime();
          algo.executeSegmentation(inSeedsSequence.getValue().getROIs(ROI2D.class));
          endTime = System.nanoTime();
          System.out.println("Time to run algorithm: " + (endTime - startTime)/1000000 + "ms");
          
          MessageDialog.showDialog(String.format("Power Watershed Plugin is working fine and using %s!", algo));
          
          // Show results
          Sequence result = algo.getSegmentation();
          addSequence(result);
        } catch (BadHistogramParameters e) {
          MessageDialog.showDialog(String.format("Power Watershed Plugin failed for the algorithm \"%s\"!", algo));
          e.printStackTrace();
        }
		
		
		
		
	}

	/* (non-Javadoc)
	 * @see plugins.adufour.ezplug.EzPlug#clean()
	 */
	@Override
	public void clean() {
		System.gc();
	}
	
	/**
	 * Validates the input data of the plugin.
	 * @return 0 if input data is correct, else an integer different than 0.
	 */
	private int validateInput() {
		if (inSequence.getValue() == null) {
			MessageDialog.showDialog("Error", 
					"Please select a sequence before starting the algorithm", 
					MessageDialog.ERROR_MESSAGE);
			return 1;
		}
		
		if (inSeedsSequence.isEnabled()) {
			if (inSeedsSequence.getValue() == null) {
				MessageDialog.showDialog("Error", 
						"Please select a seed sequence before starting the algorithm", 
						MessageDialog.ERROR_MESSAGE);
				return 2;
			} else
				if (inSeedsSequence.getValue().getROIs().isEmpty()) {
					MessageDialog.showDialog("Error", 
							"Please select a sequence with ROI's specifying the seeds " +
							"before starting the algorithm", 
							MessageDialog.ERROR_MESSAGE);
					return 3;
				} else {
					List<ROI> rois = inSeedsSequence.getValue().getROIs(ROI2D.class);
					if (rois.size() != 2) {
						MessageDialog.showDialog("Error", 
								"Please select a sequence with ROI's specifying the seeds " +
								"with two classes before starting the algorithm", 
								MessageDialog.ERROR_MESSAGE);
						return 4;
					}
				}
		}
		
		return 0;
  }
}
