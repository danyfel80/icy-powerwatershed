/**
 * 
 */
package algorithms.danyfel80.segmentation.graphcut;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import algorithms.danyfel80.segmentation.SegmentationAlgorithm;
import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceDataIterator;
import icy.sequence.SequenceUtil;
import icy.type.DataIteratorUtil;
import icy.type.DataType;
import icy.type.TypeUtil;
import plugins.kernel.roi.roi3d.ROI3DArea;
import plugins.ylemontag.histogram.BadHistogramParameters;
import plugins.ylemontag.histogram.Histogram;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
// TODO fix for color images
public class GraphCutSegmentation extends SegmentationAlgorithm{

  private Sequence inSequence;
  private Sequence treatedSequence;
  private List<ROI> seeds;
  private double lambda;
  private boolean use8Connected;

  private int sx, sy, sz;
  private int numE, numN;
  private GraphCutMethod graphCut;

  private double edgeVariance;

  private Sequence resultSegmentation;
  private List<? extends ROI> resultROIs;

  /**
   * 
   * @param sequence
   * @param seeds
   * @param lambda
   * @param use8Connected
   * @throws BadHistogramParameters
   */
  public GraphCutSegmentation(Sequence sequence, List<ROI> seeds, double lambda, double edgeVariance, Boolean use8Connected) throws BadHistogramParameters {
    inSequence = sequence;
    this.seeds = seeds;
    this.lambda = lambda;
    this.use8Connected = use8Connected;
    this.edgeVariance = edgeVariance;
    /*if (inSequence.getDataType_() == DataType.BYTE) {
      treatedSequence = inSequence;
    } else {
      treatedSequence = SequenceUtil.convertToType(inSequence, DataType.BYTE, true);
    }*/
    treatedSequence = SequenceUtil.toGray(inSequence);
    treatedSequence.setName(inSequence.getName());
    sx = treatedSequence.getSizeX();
    sy = treatedSequence.getSizeY();
    sz = treatedSequence.getSizeZ();

    prepareGraph(seeds);
    System.out.println("Graph Cut Segmentation(lambda=" + lambda + ", variance=" + edgeVariance + ", " + (use8Connected?"8": "4") + "-connect)");
  }

  @Override
  protected void prepareGraph(List<ROI> seeds) throws BadHistogramParameters{

    int x, y, z, dx, dy, dz, currPos, neighPos, i;
    double probaSrc, probaSnk;

    // 1. Graph Instantiation
    int sxy = sx*sy;
    numN = sxy*sz;
    numE = (use8Connected)?
        sz*(12*sxy - 9*sx - 9*sy + 4) - sy*(8*sx + 6) + 6*sx - 2:
          3*sxy*sz - sxy - sx*sz - sy*sz;

    graphCut = new GraphCutMethod(numN, numE);

    // 2. Neighbor node weights specification
    byte[][] seqData = treatedSequence.getDataXYZAsByte(0, 0);
    //edgeVariance = computeVariance();
    // - For each pixel
    for (z = 0; z < sz; z++) {
      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          currPos = z*sxy + y*sx + x;

          // - Take neighbors
          for (dx = 0; dx < 2 && dx+x < sx; dx++) {
            for (dy = 0; dy < 2 && dy+y < sy; dy++) {
              for (dz = 0; dz < 2 && dz+z < sz; dz++) {
                if (dx+dy+dz > 0) {
                  if(use8Connected || dx+dy == 0|| dx+dz == 0 || dy+dz == 0) {
                    neighPos = (z+dz)*sxy + (y+dy)*sx + (x+dx);
                    
                    graphCut.setEdgeWeight(
                        currPos, 
                        neighPos, 
                        (float)getEdgeLikelihood(
                            TypeUtil.unsign(seqData[z][y*sx + x]), 
                            TypeUtil.unsign(seqData[z+dz][(y+dy)*sx + (x+dx)]), 
                            dx, dy, dz)
                        ); // Sets weight on edges in both directions 
                  }
                }
              }
            }
          }
        }
      }
    }

    // 3. Terminal node weights specification
    Histogram[] seedHistograms = getSeedHistograms();
    double[] seedSizes = new double[seeds.size()];
    for (i = 0; i < seedSizes.length; i++) {
      seedSizes[i] = seeds.get(i).getNumberOfPoints();
    }
    for (z = 0; z < sz; z++) {
      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          currPos = z*sxy + y*sx + x;

          probaSrc = -Math.log((double)(seedHistograms[0].getBin(TypeUtil.unsign(seqData[z][y*sx + x])/4).getCount()) / (double)(seedSizes[0]));
          probaSrc = (probaSrc > graphCut.maxTerminalWeight)? graphCut.maxTerminalWeight: probaSrc;
          //          System.out.print("proba " + probaSrc);
          probaSnk = -Math.log((double)(seedHistograms[1].getBin(TypeUtil.unsign(seqData[z][y*sx + x])/4).getCount()) / (double)(seedSizes[1]));
          probaSnk = (probaSnk > graphCut.maxTerminalWeight)? graphCut.maxTerminalWeight: probaSnk;
          //          System.out.println(" " + probaSnk);
          //          probaSnk = -Math.log(1.0 - (seedHistograms[0].getBin(TypeUtil.unsign(seqData[z][y*sx + x])).getCount() / (double)seedSizes[0]));
          graphCut.setTerminalWeights(currPos, (float)(lambda*probaSrc), (float)(lambda*probaSnk), false);
        }
      }
    }
    SequenceDataIterator sDI;
    for (i = 0; i < 2; i++) {
       sDI = new SequenceDataIterator(treatedSequence, seeds.get(i));
      while (!sDI.done()) {
        currPos = sDI.getPositionZ()*sxy + sDI.getPositionY()*sx + sDI.getPositionX();
        graphCut.setTerminalWeights(currPos, graphCut.maxTerminalWeight*i, graphCut.maxTerminalWeight*(1-i), false);
        sDI.next();
      }
    }

  }

  /**
   * Computes the gray level variance of the treated sequence
   * @return gray level variance 
   */
  @SuppressWarnings("unused")
  private double computeVariance() {
    int z, xy;
    
    byte[][] seqData = treatedSequence.getDataXYZAsByte(0, 0);
    int sxy =  sx*sy;
    double mean = 0;
    for (z = 0; z < sz; z++) {
      for (xy = 0; xy < sxy; xy++) {
        mean += seqData[z][xy];
      }
    }
    mean /= sxy*sz;

    double variance = 0;
    double val;
    for (z = 0; z < sz; z++) {
      for (xy = 0; xy < sxy; xy++) {
        val = mean - seqData[z][xy];
        variance += val*val;
      }
    }
    variance /= sxy*sz;
    return variance;
  }

  private double getEdgeLikelihood(int value1, int value2, int dx, int dy, int dz) {
    int diff = value2 - value1;
    return Math.exp(-(diff*diff)/(2*edgeVariance))/
        Math.sqrt(dx*dx + dy*dy + dz*dz);
  }

  private Histogram[] getSeedHistograms() throws BadHistogramParameters {
    Histogram[] hists = new Histogram[seeds.size()];
    for (int i = 0; i < seeds.size(); i++) {
      hists[i] = Histogram.compute(treatedSequence, seeds.get(i), true, 64, 0, 255);
    }
    return hists;
  }

  @Override
  public void executeSegmentation() {
    int x, y, z, r;
    int sxy = sx*sy;

    double maxFlow = graphCut.computeMaximumFlow(false, null);
    System.out.printf("Max flow = %f%n", maxFlow);

    List<ROI3DArea> resROI = new ArrayList<ROI3DArea>();
    for (r = 0; r < 2; r++) {
      resROI.add(new ROI3DArea());
    }

    for (z = 0; z < sz; z++) {

      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          resROI.get((graphCut.getTerminal(z*sxy + y*sx + x) == Terminal.BACKGROUND)? 0: 1).addPoint(x, y, z);
        }
      }
    }
    resROI.get(0).setColor(seeds.get(0).getColor());
    resROI.get(1).setColor(seeds.get(1).getColor());
    this.resultROIs = resROI;
  }

  @Override
  public Sequence getSegmentationSequence() {
    if (resultSegmentation == null) {
      resultSegmentation = new Sequence();
      resultSegmentation.beginUpdate();
      try
      {
        for (int z = 0; z < sz; z++)
          resultSegmentation.setImage(0, z, new IcyBufferedImage(sx, sy, 1, DataType.UBYTE));

        // set value from ROI(s)
        int val = 0;
        for (ROI roi : resultROIs){
          if (!roi.getBounds5D().isEmpty())
            DataIteratorUtil.set(new SequenceDataIterator(resultSegmentation, roi), val);
          val++;
        }

        // notify data changed
        resultSegmentation.dataChanged();
      }
      finally
      {
        resultSegmentation.endUpdate();
      }
      resultSegmentation.setName(inSequence.getName() + 
          "(var" + edgeVariance + 
          ", lambda" + lambda + 
          ", " + ((use8Connected)? "8": "4") + "-conn)");
    }
    return resultSegmentation;
  }

  @Override
  public Collection<? extends ROI> getSegmentationROIs() {
    return resultROIs;
  }

  @SuppressWarnings("unchecked")
  @Override
  public Sequence getSegmentationSequenceWithROIs() {
    Sequence result = SequenceUtil.getCopy(treatedSequence);
    result.addROIs((Collection<ROI>) resultROIs, false);
    result.setName(inSequence.getName() + 
        "Segmentation(var" + edgeVariance + 
        ", lambda" + lambda + 
        ", " + ((use8Connected)? "8": "4") + "-conn)");
    return result;
  }

  public Sequence getTreatedSequence() {
    return this.treatedSequence;
  }

}
