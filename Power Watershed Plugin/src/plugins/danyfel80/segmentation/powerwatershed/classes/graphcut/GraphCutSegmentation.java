package plugins.danyfel80.segmentation.powerwatershed.classes.graphcut;

import icy.roi.BooleanMask2D;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import icy.type.point.Point5D;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm;
import plugins.danyfel80.segmentation.powerwatershed.classes.graphcut.Graph.TerminalType;
import plugins.kernel.roi.roi2d.ROI2DArea;
import plugins.ylemontag.histogram.BadHistogramParameters;
import plugins.ylemontag.histogram.Histogram;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class GraphCutSegmentation extends SegmentationAlgorithm {

  private final double lambda;
  private final double K;

  @SuppressWarnings("unused")
  private Sequence inSequence;
  private Sequence inGraySequence;
  private int sizeX;
  private int sizeY;
  private int sizeZ;

  private Graph graph;

  private List<Color> colors;
  private Sequence segSequence;

  /**
   * @param lambda 
   * @param inSequence 
   * 
   */
  public GraphCutSegmentation(Sequence inSequence, double lambda) {
    this.inSequence = inSequence;
    this.inGraySequence = SequenceUtil.toGray(inSequence);
    this.inGraySequence = SequenceUtil.convertToType(inSequence, DataType.DOUBLE, false);

    this.lambda = lambda;
    this.K = 27.0;

    this.sizeX = inSequence.getSizeX();
    this.sizeY = inSequence.getSizeY();
    this.sizeZ = inSequence.getSizeZ();
    prepareGraph();

    segSequence = SequenceUtil.getCopy(inSequence);
//    segSequence = new Sequence(inSequence.getName() + "_Segmentation");
//    segSequence.beginUpdate();
//    for (int z = 0; z < sizeZ; z++) {
//      segSequence.setImage(0, z, new IcyBufferedImage(sizeX, sizeY, 3, DataType.DOUBLE));
//    }
//    segSequence.endUpdate();
  }

  /* (non-Javadoc)
   * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#prepareGraph()
   */
  @Override
  protected void prepareGraph() {
    if (sizeZ == 1)
      graph = new Graph(sizeX*sizeY, sizeX*sizeY*8);
    else
      graph = new Graph(sizeX*sizeY*sizeZ, sizeX*sizeY*sizeZ*26);

    graph.addNodes(sizeX*sizeY*sizeZ);

    double variance = calculateVariance();

    // set neighbor edges
    double[][][] inGrayData = inGraySequence.getDataXYCZAsDouble(0);
    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {

          int idP = getNodeId(x, y, z);
          for (int dz = -1; dz <= 1; dz++) {
            for (int dx = -1; dx <= 1; dx++) {
              for (int dy = -1; dy <= 1; dy++) {

                if (dz != 0 || dx!= 0 || dy != 0) {
                  if (z+dz >= 0 && z+dz < sizeZ &&
                      x+dx >= 0 && x+dx < sizeX &&
                      y+dy >= 0 && y+dy < sizeY) {
                    int idQ = getNodeId(x+dx, y+dy, z+dz);
                    if (idQ > idP) {
                      double diff = inGrayData[z][0][x + y*sizeX];
                      diff -= inGrayData[(z+dz)][0][(x+dx) + (y+dy)*sizeX]; 
                      double weight = Math.exp(-(diff*diff)/(2.0 * variance)) * (1.0 / Math.sqrt((double)(dx*dx + dy*dy + dz*dz)));
                      graph.addEdge(idP, idQ, weight, weight);
                    }
                  }
                }

              }
            }
          }

        }
      }
    }
  }

  private double calculateVariance() {
    double[][][] inGrayData = inGraySequence.getDataXYCZAsDouble(0);
    double variance = 0;
    double mean = 0;
    int amount = sizeX*sizeY*sizeZ;

    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          mean += inGrayData[z][0][x + y*sizeX];
        }
      }
    }
    mean /= amount;

    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          double val = inGrayData[z][0][x + y*sizeX] - mean;
          variance += val*val;
        }
      }
    }
    variance /= amount;

    return variance;
  }

  private int getNodeId(int x, int y, int z) {
    return z*sizeX*sizeY + y*sizeX + x;
  }



  /* (non-Javadoc)
   * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#executeSegmentation(java.util.List)
   */
  @Override
  public void executeSegmentation(List<ROI> seeds) throws BadHistogramParameters {
    // calculate terminal edges
    colors = new ArrayList<>();
    for (ROI roi : seeds) {
      colors.add(roi.getColor());
      System.out.println(roi.getColor());
    }
    //paintSeeds(seeds);
    addTerminalEdges(seeds);

    double maxFlow = graph.computeMaxFlow();
    System.out.println(maxFlow);

    paintSegmentation();
  }

  //  private void paintSeeds(List<ROI> seeds) {
  //    segSequence.beginUpdate();
  //    double[][][] segData = segSequence.getDataXYCZAsDouble(0);
  //    for (int z = 0; z < sizeZ; z++) {
  //      for (int x = 0; x < sizeX; x++) {
  //        for (int y = 0; y < sizeY; y++) {
  //          if (seeds.get(0).contains(new Point5D.Double(x,y,z,0,0)))
  //            segData[z][0][x*sizeY + y] = 1;
  //          else if (seeds.get(1).contains(new Point5D.Double(x,y,z,0,0)))
  //            segData[z][0][x*sizeY + y] = 2;
  //        }
  //      }
  //    }
  //    segSequence.dataChanged();
  //    segSequence.endUpdate();
  //  }

  /**
   * Computes the histogram for each of the seed ROIs.
   * @return Histogram array (one histogram for each ROI).
   * @throws BadHistogramParameters 
   */
  private Histogram[] getClassesHistograms(List<ROI> seeds) throws BadHistogramParameters {
    Histogram[] histograms = new Histogram[seeds.size()];
    for (int i = 0; i < seeds.size(); i++) {
      histograms[i] = Histogram.compute(inGraySequence, seeds.get(i), true, 256, 0, 255);
    }
    return histograms;
    /*
    int [] amount = new int[2];
    double[][] histograms = new double[2][256];

    double [][][] segData = segSequence.getDataXYCZAsDouble(0); 
    double [][][] imageData = inGraySequence.getDataXYCZAsDouble(0);

    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          double v = segData[z][0][x + y*sizeX];
          if (v > 0) {
            int p = ((int)imageData[z][0][x + y*sizeX]);
            histograms[(int)v-1][p] += 1.0;
            amount[(int)v-1]++;
          }
        }
      }
    }
    double maxVal=0, sum0=0, sum1=0;

    for (int i = 0; i < 256; i++) {
      //histograms[0][i] /= (double)amount[0];    
      //histograms[1][i] /= (double)amount[1];
      sum0 += histograms[0][i];
      sum1 += histograms[1][i];
      maxVal = Math.max(histograms[0][i], histograms[1][i]);
      System.out.printf("[%d]=%f\n", i, histograms[0][i]);
    }
    System.out.println("Max histo val = " + maxVal);
    System.out.println("sum histo 0 = " + sum0);
    System.out.println("sum histo 1 = " + sum1);
    System.out.println("quant histo 0 = " + amount[0]);
    System.out.println("quant histo 1 = " + amount[1]);
    return histograms;
     */
  }

  private void addTerminalEdges(List<ROI> seeds) throws BadHistogramParameters {
    double [][][] grayData = inGraySequence.getDataXYCZAsDouble(0);
    Histogram[] histograms = getClassesHistograms(seeds);
    double [] roiSizes = new double[seeds.size()];

    for (int i = 0; i < seeds.size(); i++) {
      roiSizes[i] = seeds.get(i).getNumberOfPoints();
    }

    //    for (int i = 0; i < 256; i++) {
    //      System.out.println("[" + i + "] = " + (histograms[0].getBin(i).getCount() / roiSizes[0]));
    //    }

    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          double capacitySource = 0.0;
          double capacitySink = 0.0;
          if (seeds.get(0).contains(new Point5D.Double(x,y,z,0,0))) { // source seed
            capacitySource = 0;
            capacitySink = K;
          } else if (seeds.get(0).contains(new Point5D.Double(x,y,z,0,0))) { // sink seed
            capacitySource = K;
            capacitySink = 0;
          } else { // if not seed, then use histograms
            int p = (int)grayData[z][0][x + y*sizeX];
            double prPSrc = histograms[0].getBin(p).getCount() / roiSizes[0];
            double prPSnk = histograms[1].getBin(p).getCount() / roiSizes[1];
            capacitySource = lambda * ((prPSrc > 0.000001)? (-Math.log(prPSrc)): K);
            capacitySink = lambda * ((prPSnk > 0.000001)? (-Math.log(prPSnk)): K);
          }

          graph.addTerminalWeights(getNodeId(x, y, z), capacitySink, capacitySource);
        }
      }
    }
  }

  private void paintSegmentation() {
    //segSequence.beginUpdate();
    //double[][][] segData = segSequence.getDataXYCZAsDouble(0);
    
    
    for (int z = 0; z < sizeZ; z++) {
      
      boolean[] maskSource = new boolean[sizeX*sizeY];
      boolean[] maskSink = new boolean[sizeX*sizeY];
      
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          if (graph.getSegment(getNodeId(x, y, z)) == TerminalType.SOURCE) {
            maskSource[x + y*sizeX] = true;
            maskSink[x + y*sizeX] = false;
          }
          if (graph.getSegment(getNodeId(x, y, z)) == TerminalType.SINK) {
            maskSource[x + y*sizeX] = false;
            maskSink[x + y*sizeX] = true;
          }

          //          if (color < 2) {
          //            System.out.println(color);
          //            segData[z][0][x + y*sizeX] = colors.get(color).getRed();
          //            segData[z][1][x + y*sizeX] = colors.get(color).getGreen();
          //            segData[z][2][x + y*sizeX] = colors.get(color).getBlue();
          //          }
          //          else {
//          segData[z][0][x + y*sizeX] = 1;
//          segData[z][1][x + y*sizeX] = 1;
//          segData[z][2][x + y*sizeX] = 1;
          //          }
        }
      }
      
      BooleanMask2D mask2Dsource = new BooleanMask2D(segSequence.getBounds2D(), maskSource);
      ROI2DArea roiSourceZ = new ROI2DArea(mask2Dsource);
      roiSourceZ.setZ(z);
      roiSourceZ.setColor(colors.get(0));
      segSequence.addROI(roiSourceZ);
      BooleanMask2D mask2Dsink = new BooleanMask2D(segSequence.getBounds2D(), maskSink);
      ROI2DArea roiSinkZ = new ROI2DArea(mask2Dsink);
      roiSinkZ.setZ(z);
      roiSinkZ.setColor(colors.get(1));
      segSequence.addROI(roiSinkZ);
    }

    //segSequence.dataChanged();
    //segSequence.endUpdate();
  }

  /* (non-Javadoc)
   * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#getSegmentation()
   */
  @Override
  public Sequence getSegmentation() {
    return segSequence;
  }

}
