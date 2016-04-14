package plugins.danyfel80.segmentation.powerwatershed;

import java.util.ArrayList;
import java.util.List;

import algorithms.danyfel80.segmentation.SegmentationAlgorithm;
import algorithms.danyfel80.segmentation.graphcut.GraphCutSegmentation;
import algorithms.danyfel80.segmentation.powerwatershed.PowerWatershedSegmentation;
import icy.gui.dialog.MessageDialog;
import icy.roi.ROI;
import icy.sequence.Sequence;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVar;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarEnum;
import plugins.adufour.ezplug.EzVarListener;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.ylemontag.histogram.BadHistogramParameters;

/**
 * @author Daniel Felipe Gonzalez Obando
 * 
 *         Graph Segmentation Framework Plugin
 * 
 *         This plugin implements the segmentation framework proposed in the
 *         article
 *         "Power Watershed: A Unifying Graph-Based Optimization Framework" by
 *         Camille Couprie, Leo Grady, Laurent Najman, Hugues Talbot.
 */
public class GraphSegmentation extends EzPlug {

  // General input
  private EzVarSequence inSequence;
  private Sequence inSeq;
  private EzVarSequence inSeedsSequence;
  private Sequence inSeeds;
  private EzVarEnum<SegmentationType> inSegmentationType;
  
  // Specific input
  // - Graphcuts
  private EzVarBoolean inUse8Connected;
  private EzVarDouble inLambda;
  private EzVarDouble inEdgeVariance;
  private EzVarBoolean inShowProbas;
  
  // - PowerWatershed Q2
  private EzVarBoolean inUseGrayLevels;
  private EzVarBoolean inUseGeodesicReconstruction;

  /* (non-Javadoc)
   * @see plugins.adufour.ezplug.EzPlug#initialize()
   */
  @Override
  protected void initialize() {

    // Input instantiation and setup
    inSequence = new EzVarSequence("Sequence");
    inSequence.setToolTipText("The sequence to be segmented.");
    inSequence.addVarChangeListener(new EzVarListener<Sequence>() {
      
      @Override
      public void variableChanged(EzVar<Sequence> source, Sequence newValue) {
        if(!inSeedsSequence.isEnabled()) {
          inSeedsSequence.setValue(inSequence.getValue());
        }
      }
    });

    inSeedsSequence = new EzVarSequence("Seeds");
    inSeedsSequence.setOptional(true);
    inSeedsSequence.setToolTipText("The seeds for the segmentation, "
        + "at least two classes must be specified.");
    
    inSegmentationType = new EzVarEnum<>("Segmentation Type", SegmentationType.values());
    
    // specific params
    @SuppressWarnings("rawtypes")
    final List<EzVar> optional = new ArrayList<>();
    
    inUse8Connected = new EzVarBoolean("Use 8 connected image", false);
    inUse8Connected.setToolTipText("If true uses 8 neighbors instead of 4.");
    inLambda = new EzVarDouble("Lambda");
    inLambda.setValue(1.0);
    inLambda.setToolTipText("The importance of inter-node weight");
    inLambda.setMinValue(0.0);
    inEdgeVariance = new EzVarDouble("Edge variance");
    inEdgeVariance.setValue(25.0);
    inEdgeVariance.setToolTipText("The variance in pixels' neighborhood");
    inEdgeVariance.setMinValue(0.001);
    inEdgeVariance.setMaxValue(65025.0);
    inShowProbas = new EzVarBoolean("Show probabilities", false);
    inShowProbas.setToolTipText("If checked the gradient and terminal probabilities sequences are created.");;
    optional.add(inUse8Connected);
    optional.add(inLambda);
    optional.add(inEdgeVariance);
    optional.add(inShowProbas);
    
    inUseGrayLevels = new EzVarBoolean("Use gray levels", false);
    inUseGrayLevels.setToolTipText("If true sequence is converted to gray levels.");
    inUseGeodesicReconstruction = new EzVarBoolean("Use geodesic reconstruction.", false);
    inUseGeodesicReconstruction.setToolTipText("If true the weight is normalized using a geodesic reconstruction.");
    optional.add(inUseGrayLevels);
    optional.add(inUseGeodesicReconstruction);
    
    
    
    inSegmentationType.addVarChangeListener(new EzVarListener<SegmentationType>() {
      @Override
      public void variableChanged(EzVar<SegmentationType> source, SegmentationType newValue) {
        for (@SuppressWarnings("rawtypes") EzVar option : optional) {
          option.setVisible(false);
        }
        
        switch (newValue) {
        case CollapseToSeeds:
          break;
        case GraphCuts:
          inUse8Connected.setVisible(true);
          inLambda.setVisible(true);
          inEdgeVariance.setVisible(true);
          inShowProbas.setVisible(true);
          break;
        case PowerWatershedQ1:
          break;
        case PowerWatershedQ2:
          inUseGrayLevels.setVisible(true);
          inUseGeodesicReconstruction.setVisible(true);
          break;
        case RandomWalker:
          break;
        case ShortestPathForest:
          break;
        case l1NormVoronoi:
          break;
        case l2NormVoronoi:
          break;
        default:
          break;
        }
      }
    });

    
    

    // Input components addition to UI
    addEzComponent(inSequence);
    addEzComponent(inSeedsSequence);
    inSeedsSequence.setEnabled(false);
    addEzComponent(inSegmentationType);
    addEzComponent(inUse8Connected);
    addEzComponent(inLambda);
    addEzComponent(inEdgeVariance);
    addEzComponent(inShowProbas);
    addEzComponent(inUseGrayLevels);
    addEzComponent(inUseGeodesicReconstruction);
    
    inSegmentationType.setValue(SegmentationType.PowerWatershedQ2);
  }

  /* (non-Javadoc)
   * @see plugins.adufour.ezplug.EzPlug#execute()
   */
  @Override
  protected void execute() {

    // Get the used algorithm
    SegmentationType algoType = inSegmentationType.getValue();
    SegmentationAlgorithm algo = null;
    this.inSeq = inSequence.getValue();
    this.inSeeds = inSeedsSequence.getValue();
    
    if (!inSeedsSequence.isEnabled()) {
      inSeeds = inSeq;
    }
    
    // Input validation
    if (validateInput(algoType) != 0) {
      return;
    }

    long startTime = 0, endTime = 0;

    try {
      switch (algoType) {
      case CollapseToSeeds:

        break;
      case l2NormVoronoi:

        break;
      case l1NormVoronoi:

        break;
      case GraphCuts:
        startTime = System.nanoTime();
        GraphCutSegmentation a3 = new GraphCutSegmentation(
            inSequence.getValue(), 
            inSeedsSequence.getValue().getROIs(ROI.class),
            inLambda.getValue().floatValue(),
            inEdgeVariance.getValue(),
            inUse8Connected.getValue(),
            inShowProbas.getValue());
        endTime = System.nanoTime();
//        addSequence(a3.getTreatedSequence());
        if (inShowProbas.getValue()) {
          addSequence(a3.getGradientSequence());
          addSequence(a3.getTerminalProbabilitiesSequence());
        }
        algo = a3;
//        return;
        break;
      case RandomWalker:

        break;
      case PowerWatershedQ1:

        break;
      case PowerWatershedQ2:
        startTime = System.nanoTime();
        PowerWatershedSegmentation a6 = new PowerWatershedSegmentation(
            inSequence.getValue(), 
            inSeedsSequence.getValue().getROIs(ROI.class), 
            inUseGeodesicReconstruction.getValue(),
            inUseGrayLevels.getValue());
        endTime = System.nanoTime();
        algo = a6;
//        addSequence(a6.getTreatedSequence());
        addSequence(a6.getSequenceGradient());
//        return;
        break;
      case ShortestPathForest:

        break;
      }
    } catch (BadHistogramParameters e) {
      MessageDialog.showDialog(String.format("Power Watershed Plugin failed for the algorithm \"%s\"!", algo));
      e.printStackTrace();
    }

    System.out.println("Time to prepare algorithm: " + (endTime - startTime)/1000000 + "ms");

    // Execute algorithm

    startTime = System.nanoTime();
    algo.executeSegmentation();
    endTime = System.nanoTime();
    System.out.println("Time to run algorithm: " + (endTime - startTime)/1000000 + "ms");

    MessageDialog.showDialog(String.format("Power Watershed Plugin is working fine and using %s!", algo));

    // Show results
    Sequence result = algo.getSegmentationSequenceWithROIs();
    addSequence(algo.getSegmentationSequence());
    addSequence(result);
    System.gc();
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
  private int validateInput(SegmentationType algoType) {
    if (inSeq == null) {
      MessageDialog.showDialog("Error", 
          "Please select a sequence before starting the algorithm", 
          MessageDialog.ERROR_MESSAGE);
      return 1;
    }

    
    if (inSeeds == null) {
      MessageDialog.showDialog("Error", 
          "Please select a seed sequence before starting the algorithm", 
          MessageDialog.ERROR_MESSAGE);
      return 2;
    } else {
      if (inSeeds.getROIs().isEmpty()) {
        MessageDialog.showDialog("Error", 
            "Please select a sequence with ROI's specifying the seeds " +
                "before starting the algorithm", 
                MessageDialog.ERROR_MESSAGE);
        return 3;
      } else {
        List<ROI> rois = inSeeds.getROIs(ROI.class);
        
        if (algoType.equals(SegmentationType.GraphCuts) && rois.size() != 2) {
          MessageDialog.showDialog("Error", 
              "Please select a sequence with ROI's specifying the seeds " +
                  "with two classes before starting the algorithm", 
                  MessageDialog.ERROR_MESSAGE);
          return 4;
        } else if (rois.size() < 2) {
          MessageDialog.showDialog("Error", 
              "Please select a sequence with ROI's specifying the seeds " +
                  "with at least two classes before starting the algorithm", 
                  MessageDialog.ERROR_MESSAGE);
          return 5;
        }
      }
    }
    
    return 0;
  }
}
