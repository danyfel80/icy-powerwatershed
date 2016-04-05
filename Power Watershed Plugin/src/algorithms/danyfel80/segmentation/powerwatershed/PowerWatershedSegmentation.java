package algorithms.danyfel80.segmentation.powerwatershed;

import java.util.List;

import algorithms.danyfel80.segmentation.SegmentationAlgorithm;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import plugins.kernel.roi.roi3d.ROI3DArea;
import plugins.ylemontag.histogram.BadHistogramParameters;

/**
 * Power Watershed implementation
 * @author Daniel Felipe Gonzalez Obando
 */
public class PowerWatershedSegmentation extends SegmentationAlgorithm {

  private Sequence inSequence;
  private Sequence inGraySequence;
  private boolean useGeo;
  
  GraphPW graph;
  
  public PowerWatershedSegmentation(Sequence sequence, List<ROI> seeds, boolean useGeo) throws BadHistogramParameters {
    inSequence = sequence;
    inGraySequence = SequenceUtil.toGray(sequence);
    this.useGeo = useGeo;
    prepareGraph(seeds);
  }
  
  /* (non-Javadoc)
   * @see algorithms.danyfel80.segmentation.SegmentationAlgorithm#prepareGraph()
   */
  @Override
  protected void prepareGraph(List<ROI> seeds) throws BadHistogramParameters {
    long startTime = System.nanoTime();
    graph = new GraphPW(inGraySequence, seeds);
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
    graph.executeWatershed(useGeo);
  }

  /* (non-Javadoc)
   * @see algorithms.danyfel80.segmentation.SegmentationAlgorithm#getSegmentation()
   */
  @Override
  public Sequence getSegmentation() {
    Sequence copy = SequenceUtil.getCopy(inSequence);
    
    
    copy.beginUpdate();
    for (ROI3DArea roi: graph.getSegmentationROIs()) {
      copy.addROI(roi);
    }
    copy.endUpdate();
    return copy;
//    return graph.getSegmentation();
  }

}
