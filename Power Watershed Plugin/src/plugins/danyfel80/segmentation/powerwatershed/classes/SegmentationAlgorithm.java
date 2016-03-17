/**
 * 
 */
package plugins.danyfel80.segmentation.powerwatershed.classes;

import icy.roi.ROI;
import icy.sequence.Sequence;

import java.util.List;

/**
 * @author Daniel Felipe Gonzalez Obando
 * This abstract class represent the generalization of a segmentation algorithm.
 * In general, as all the algorithms used in this project use graphs, the
 * graph construction methods are stored in this abstract class.
 */
public abstract class SegmentationAlgorithm {
	
	
	/**
	 * Prepares the graph for the algorithm
	 * @param image The image to segment.
	 */
	protected abstract void prepareGraph();

	/**
	 * Executes the segmentation algorithm.
	 * @param list The seeds for the segmentation. At least two classes must be
	 * given.
	 */
	public abstract void executeSegmentation(List<ROI> seeds);

	public abstract Sequence getSegmentation();
	
}
