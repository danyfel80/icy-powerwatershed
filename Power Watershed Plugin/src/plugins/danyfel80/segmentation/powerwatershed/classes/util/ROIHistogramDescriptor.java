/**
 * 
 */
package plugins.danyfel80.segmentation.powerwatershed.classes.util;

import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.roi.ROIDescriptor;
import icy.sequence.Sequence;

import java.util.ArrayList;

/**
 * @author Daniel Felipe Gonzalez Obando
 * Histogram ROI descriptor class (see {@link ROIDescriptor}).
 * Gets the histogram for single-channel sequences.
 */
public class ROIHistogramDescriptor extends ROIDescriptor {

	public static final String ID = "Histogram";
	
	private int nBins;
	private int c;
	private boolean normalized;
	
	/**
	 * Default constructor.
	 * @param nBins Number of bins for the histogram
	 * @param c Channel used to create histogram. -1 uses mean value of all channels
	 */
	public ROIHistogramDescriptor(int nBins, boolean normalized, int c)
  {
      super(ID, "Histogram", ArrayList.class);
      this.nBins = nBins;
      this.normalized = normalized;
      
      if (c > -1)
      	this.c = c;
      else
      	this.c = -1;
  }

	@Override
  public String getUnit(Sequence sequence)
  {
		if (normalized) {
	    return "%";
    } else {
      return "px";
    }
  }
	
	/* (non-Javadoc)
	 * @see icy.roi.ROIDescriptor#getDescription()
	 */
	@Override
	public String getDescription() {
		return "Histogram of the ROI";
	}

	/* (non-Javadoc)
	 * @see icy.roi.ROIDescriptor#compute(icy.roi.ROI, icy.sequence.Sequence)
	 */
	@Override
	public Object compute(ROI roi, Sequence sequence)
	    throws UnsupportedOperationException {
		
		double[] histogram = new double[nBins];
		int amount = 0;
		
		int sizeX = sequence.getSizeX();
		int sizeY = sequence.getSizeY();
		int sizeZ = sequence.getSizeZ();
		int sizeC = sequence.getSizeC();
		
		double minVal;
		double maxVal;
		
		if (c != -1) {
			minVal = sequence.getChannelMin(c);
			maxVal = sequence.getChannelMax(c);
		} else {
			minVal = sequence.getChannelMin(0);
			maxVal = sequence.getChannelMax(0);
			for (int channel = 1; channel < histogram.length; channel++) {
				minVal = Math.min(sequence.getChannelMin(channel-1), sequence.getChannelMin(channel));
				maxVal = Math.max(sequence.getChannelMin(channel-1), sequence.getChannelMin(channel));
      }
		}
		
		double[][][] image = sequence.getDataXYCZAsDouble(0);
		
		
		for (int x = 0; x < sizeX; x++) {
			for (int y = 0; y < sizeY; y++) {
				for (int z = 0; z < sizeZ; z++) {
					BooleanMask2D mask = roi.getBooleanMask2D(z, 0, (c==-1? 0: c), true);
					if (mask.contains(x, y)) {
						double val = 0;
				    if (c != -1) {
				    	val = image[z][c][x + y * sizeX];
				    } else {
				    	for (int ch = 0; ch < sizeC; ch++) {
				    		val += image[z][ch][x + y * sizeX];
				    	}
				    	val /= 3;
				    }
				    
				    val -= minVal;
			    	val /= maxVal;
			    	val *= nBins;
			    	histogram[(int)val] += 1;
			    	amount += 1;
					}
		    }
	    }
    }
		
		if (normalized) {
			for (int b = 0; b < nBins; b++) {
				histogram[b] /= amount; 
			}
		}
		return histogram;
	}

}
