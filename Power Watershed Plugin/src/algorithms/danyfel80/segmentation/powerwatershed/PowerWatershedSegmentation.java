package algorithms.danyfel80.segmentation.powerwatershed;

import java.util.Collection;
import java.util.List;

import algorithms.danyfel80.segmentation.SegmentationAlgorithm;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import plugins.kernel.roi.roi3d.ROI3DArea;
import plugins.ylemontag.histogram.BadHistogramParameters;

/**
 * Power Watershed implementation
 * @author Daniel Felipe Gonzalez Obando
 */
public class PowerWatershedSegmentation extends SegmentationAlgorithm {

  private Sequence inSequence;
  private Sequence treatedSequence;
  private boolean useGeo;
  private boolean isGraySequence;
  
  GraphPW graph;
  
  public PowerWatershedSegmentation(Sequence sequence, List<ROI> seeds, boolean useGeo, boolean useGrayLevels) throws BadHistogramParameters {
    this.useGeo = useGeo;
    inSequence = sequence;
    double max  = 0;
    for (int i = 0; i < inSequence.getSizeC(); i++) {
      max = Math.max(inSequence.getChannelMax(0), max);
    }
    treatedSequence = (max > 255)? SequenceUtil.convertToType(inSequence, DataType.UBYTE, true): inSequence; 
    isGraySequence = (inSequence.getSizeC() == 1)? true: useGrayLevels;
    
    if(isGraySequence) {
      treatedSequence = SequenceUtil.toGray(treatedSequence);
    }
    treatedSequence.setName(inSequence.getName());
    prepareGraph(seeds);
  }
  
  /* (non-Javadoc)
   * @see algorithms.danyfel80.segmentation.SegmentationAlgorithm#prepareGraph()
   */
  @Override
  protected void prepareGraph(List<ROI> seeds) throws BadHistogramParameters {
    long startTime = System.nanoTime();
    graph = new GraphPW(treatedSequence, seeds, useGeo, isGraySequence, false);
    long endTime = System.nanoTime();
    System.out.println("graph object created: " + ((endTime-startTime)/1000000) + " msec...");
    startTime = System.nanoTime();
    graph.SetupSeeds();
    endTime = System.nanoTime();
    System.out.println("seeds prepared: " + ((endTime-startTime)/1000000) + " msec...");
    startTime = System.nanoTime();
    graph.SetupEdges();
    endTime = System.nanoTime();
    System.out.println("edges prepared: " + ((endTime-startTime)/1000000) + " msec...");
    startTime = System.nanoTime();
    graph.calculateWeights();
    endTime = System.nanoTime();
    System.out.println("weights computed: " + ((endTime-startTime)/1000000) + " msec...");
  }

  /* (non-Javadoc)
   * @see algorithms.danyfel80.segmentation.SegmentationAlgorithm#executeSegmentation(java.util.List)
   */
  @Override
  public void executeSegmentation(){
    graph.executeWatershed();
  }

  @Override
  public Sequence getSegmentationSequence() {
    return graph.getSegmentation();
  }

  @Override
  public Collection<? extends ROI> getSegmentationROIs() {
    return (Collection<? extends ROI>)graph.getSegmentationROIs();
  }

  @Override
  public Sequence getSegmentationSequenceWithROIs() {
    Sequence copy = SequenceUtil.getCopy(treatedSequence);
    copy.setName(inSequence.getName());
    
    copy.beginUpdate();
    for (ROI3DArea roi: graph.getSegmentationROIs()) {
      copy.addROI(roi);
    }
    copy.endUpdate();
    
    return copy;
  }
  
  

}
