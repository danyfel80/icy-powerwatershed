package algorithms.danyfel80.segmentation;

import icy.roi.ROI;
import icy.sequence.Sequence;

import java.util.Collection;
import java.util.List;

import plugins.ylemontag.histogram.BadHistogramParameters;

/**
 * This abstract class represent the generalization of a segmentation algorithm.
 * In general, as all the algorithms used in this project use graphs, the graph
 * construction methods are stored in this abstract class.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public abstract class SegmentationAlgorithm {

  /**
   * Prepares the graph for the algorithm.
   * 
   * @param list
   *          The seeds for the segmentation. At least two classes must be
   *          given.
   * @throws BadHistogramParameters
   */
  protected abstract void prepareGraph(List<ROI> seeds)
      throws BadHistogramParameters;

  /**
   * Executes the segmentation algorithm.
   */
  public abstract void executeSegmentation();

  /**
   * @return Sequence with segments in colors.
   */
  public abstract Sequence getSegmentationSequence();

  /**
   * @return collection of ROIs with the found segments.
   */
  public abstract Collection<? extends ROI> getSegmentationROIs();

  /**
   * @return Sequence with segments as ROIs.
   */
  public abstract Sequence getSegmentationSequenceWithROIs();

}
