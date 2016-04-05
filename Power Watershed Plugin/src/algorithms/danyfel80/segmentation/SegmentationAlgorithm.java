/**
 * 
 */
package algorithms.danyfel80.segmentation;

import icy.roi.ROI;
import icy.sequence.Sequence;

import java.util.List;

import plugins.ylemontag.histogram.BadHistogramParameters;

/**
 * @author Daniel Felipe Gonzalez Obando
 * This abstract class represent the generalization of a segmentation algorithm.
 * In general, as all the algorithms used in this project use graphs, the
 * graph construction methods are stored in this abstract class.
 */
public abstract class SegmentationAlgorithm {
	
	
	/**
	 * Prepares the graph for the algorithm.
	 * @param list The seeds for the segmentation. At least two classes must be
     * given.
     * @throws BadHistogramParameters 
	 */
	protected abstract void prepareGraph(List<ROI> seeds) throws BadHistogramParameters;

	/**
	 * Executes the segmentation algorithm.
	 */
	public abstract void executeSegmentation();

	public abstract Sequence getSegmentation();
	
}
