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
public class GraphCutSegmentation extends SegmentationAlgorithm {

  private Sequence            inSequence;
  private Sequence            treatedSequence;
  private Sequence            graySequence;
  private List<ROI>           seeds;
  private float               lambda;
  private boolean             use8Connected;
  private boolean             createProbas;

  private int                 sx, sy, sz, sc;
  private int                 numE, numN;
  private GraphCutMethod      graphCut;

  private Histogram[]         seedHistograms;
  private double[]            seedSizes;

  private double              edgeVariance;
  // private int edgeStdDev;

  private Sequence            gradientSequence;
  private Sequence            terminalSequence;
  private Sequence            resultSegmentation;
  private List<? extends ROI> resultROIs;

  /**
   * 
   * @param sequence
   * @param seeds
   * @param lambda
   * @param use8Connected
   * @throws BadHistogramParameters
   */
  public GraphCutSegmentation(Sequence sequence, List<ROI> seeds, float lambda,
      double edgeVariance, boolean use8Connected, boolean createProbas)
      throws BadHistogramParameters {
    inSequence = sequence;
    this.seeds = seeds;
    this.lambda = lambda;
    this.use8Connected = use8Connected;
    // this.edgeVariance = edgeVariance;
    // this.edgeStdDev = (int) Math.ceil(Math.sqrt(edgeVariance));
    this.createProbas = createProbas;

    if (inSequence.getDataType_() == DataType.UBYTE) {
      treatedSequence = inSequence;
    } else {
      treatedSequence =
          SequenceUtil.convertToType(inSequence, DataType.UBYTE, true);
    }
    treatedSequence.setName(inSequence.getName());
    graySequence = SequenceUtil.toGray(inSequence);

    sx = treatedSequence.getSizeX();
    sy = treatedSequence.getSizeY();
    sz = treatedSequence.getSizeZ();
    sc = treatedSequence.getSizeC();

    this.edgeVariance = (edgeVariance != -1)? edgeVariance: computeVariance();
    System.out.println("Graph Cut Segmentation(lambda=" + lambda + ", variance="
        + this.edgeVariance + ", " + (use8Connected? "8": "4") + "-connect)");
    prepareGraph(seeds);

  }

  @Override
  protected void prepareGraph(List<ROI> seeds) throws BadHistogramParameters {

    int x, y, z, c, dx, dy, dz, currPos, neighPos, neighPos1, i, numNeigh;
    int[] currVals = new int[sc], neighVals = new int[sc],
        neighVals1 = new int[sc];
    float probaSrc, probaSnk;

    if (createProbas) {
      gradientSequence = new Sequence(inSequence.getName() + "Gradient");
      terminalSequence =
          new Sequence(inSequence.getName() + "TerminalProbabilities");
    }
    // 1. Graph Instantiation
    int sxy = sx * sy;
    numN = sxy * sz;
    numE = (use8Connected)? 13 * numN - 9 * sxy - 9 * sx * sz - 9 * sy * sz
        + 6 * sx + 6 * sy + 6 * sz - 4: 3 * sxy * sz - sxy - sx * sz - sy * sz;

    graphCut = new GraphCutMethod(numN, numE);

    // 2. Neighbor node weights specification
    float maxDiff = 0, diffN, weight;
    if (createProbas)
      gradientSequence.beginUpdate();
    byte[][][] seqData = treatedSequence.getDataXYCZAsByte(0);
    // edgeVariance = computeVariance();
    // - For each pixel
    for (z = 0; z < sz; z++) {
      float[] gradientData = null;
      if (createProbas) {
        gradientSequence.setImage(0, z,
            new IcyBufferedImage(sx, sy, 1, DataType.FLOAT));
        gradientData = gradientSequence.getDataXYAsFloat(0, z, 0);
      }
      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          currPos = z * sxy + y * sx + x;
          for (c = 0; c < sc; c++) {
            currVals[c] = TypeUtil.unsign(seqData[z][c][y * sx + x]);
          }
          diffN = 0;
          numNeigh = 0;
          // - Take neighbors
          for (dx = 0; dx < 2 && dx + x < sx; dx++) {
            for (dy = 0; dy < 2 && dy + y < sy; dy++) {
              for (dz = 0; dz < 2 && dz + z < sz; dz++) {
                if (dx + dy + dz > 0) {
                  if (use8Connected || dx + dy == 0 || dx + dz == 0
                      || dy + dz == 0) {
                    neighPos = (z + dz) * sxy + (y + dy) * sx + (x + dx);
                    for (c = 0; c < sc; c++) {
                      neighVals[c] = TypeUtil
                          .unsign(seqData[z + dz][c][(y + dy) * sx + (x + dx)]);
                    }

                    weight = (float) getEdgeLikelihood(currVals, neighVals, dx,
                        dy, dz);
                    graphCut.setEdgeWeight(currPos, neighPos, weight);
                    // Sets weight on edges in both directions
                    diffN += weight;
                    numNeigh++;
                    if (createProbas)
                      gradientData[y * sx + x] += weight;
                  }
                  if (use8Connected) {
                    if (dz == 0 && dx + dy == 2) {
                      neighPos = z * sxy + y * sx + (x + dx);
                      neighPos1 = z * sxy + (y + dy) * sx + x;
                      for (c = 0; c < sc; c++) {
                        neighVals[c] =
                            TypeUtil.unsign(seqData[z][c][y * sx + (x + dx)]);
                        neighVals1[c] =
                            TypeUtil.unsign(seqData[z][c][(y + dy) * sx + x]);
                      }
                      weight = (float) getEdgeLikelihood(neighVals, neighVals1,
                          dx, dy, dz);
                      graphCut.setEdgeWeight(neighPos, neighPos1, weight);
                      // Sets weight on edges in both directions
                      diffN += weight;
                      numNeigh++;
                      if (createProbas)
                        gradientData[y * sx + x] += weight;
                    }
                    if (dy == 0 && dx + dz == 2) {
                      neighPos = z * sxy + y * sx + (x + dx);
                      neighPos1 = (z + dz) * sxy + y * sx + x;
                      for (c = 0; c < sc; c++) {
                        neighVals[c] =
                            TypeUtil.unsign(seqData[z][c][y * sx + (x + dx)]);
                        neighVals1[c] =
                            TypeUtil.unsign(seqData[z + dz][c][y * sx + x]);
                      }
                      weight = (float) getEdgeLikelihood(neighVals, neighVals1,
                          dx, dy, dz);
                      graphCut.setEdgeWeight(neighPos, neighPos1, weight);
                      // Sets weight on edges in both directions
                      diffN += weight;
                      numNeigh++;
                      if (createProbas)
                        gradientData[y * sx + x] += weight;
                    }
                    if (dx == 0 && dy + dz == 2) {
                      neighPos = z * sxy + (y + dy) * sx + x;
                      neighPos1 = (z + dz) * sxy + y * sx + x;
                      for (c = 0; c < sc; c++) {
                        neighVals[c] =
                            TypeUtil.unsign(seqData[z][c][(y + dy) * sx + x]);
                        neighVals1[c] =
                            TypeUtil.unsign(seqData[z + dz][c][y * sx + x]);
                      }
                      weight = (float) getEdgeLikelihood(neighVals, neighVals1,
                          dx, dy, dz);
                      graphCut.setEdgeWeight(neighPos, neighPos1, weight);
                      // Sets weight on edges in both directions
                      diffN += weight;
                      numNeigh++;
                      if (createProbas)
                        gradientData[y * sx + x] += weight;
                    }

                    if (dx + dy + dz == 3) {
                      // xy z
                      neighPos = z * sxy + (y + dy) * sx + (x + dx);
                      neighPos1 = (z + dz) * sxy + y * sx + x;
                      for (c = 0; c < sc; c++) {
                        neighVals[c] = TypeUtil
                            .unsign(seqData[z][c][(y + dy) * sx + (x + dx)]);
                        neighVals1[c] =
                            TypeUtil.unsign(seqData[z + dz][c][y * sx + x]);
                      }
                      weight = (float) getEdgeLikelihood(neighVals, neighVals1,
                          dx, dy, dz);
                      graphCut.setEdgeWeight(neighPos, neighPos1, weight);
                      // Sets weight on edges in both directions
                      diffN += weight;
                      numNeigh++;
                      if (createProbas)
                        gradientData[y * sx + x] += weight;

                      // xz y
                      neighPos = (z + dz) * sxy + y * sx + (x + dx);
                      neighPos1 = z * sxy + (y + dy) * sx + x;
                      for (c = 0; c < sc; c++) {
                        neighVals[c] = TypeUtil
                            .unsign(seqData[z + dx][c][y * sx + (x + dx)]);
                        neighVals1[c] =
                            TypeUtil.unsign(seqData[z][c][(y + dy) * sx + x]);
                      }
                      weight = (float) getEdgeLikelihood(neighVals, neighVals1,
                          dx, dy, dz);
                      graphCut.setEdgeWeight(neighPos, neighPos1, weight);
                      // Sets weight on edges in both directions
                      diffN += weight;
                      numNeigh++;
                      if (createProbas)
                        gradientData[y * sx + x] += weight;

                      // yz x
                      neighPos = (z + dz) * sxy + (y + dy) * sx + x;
                      neighPos1 = z * sxy + y * sx + (x + dx);
                      for (c = 0; c < sc; c++) {
                        neighVals[c] = TypeUtil
                            .unsign(seqData[z + dz][c][(y + dy) * sx + x]);
                        neighVals1[c] =
                            TypeUtil.unsign(seqData[z][c][y * sx + (x + dx)]);
                      }
                      weight = (float) getEdgeLikelihood(neighVals, neighVals1,
                          dx, dy, dz);
                      graphCut.setEdgeWeight(neighPos, neighPos1, weight);
                      // Sets weight on edges in both directions
                      diffN += weight;
                      numNeigh++;
                      if (createProbas)
                        gradientData[y * sx + x] += weight;
                    }
                  }
                }
              }
            }
          }
          if (diffN + 1 > maxDiff)
            maxDiff = diffN + 1;
          if (createProbas)
            gradientData[y * sx + x] /= (float) numNeigh;
        }
      }
    }
    if (createProbas) {
      gradientSequence.dataChanged();
      gradientSequence.endUpdate();
    }

    // 3. Terminal node weights specification
    byte[][] grayData = graySequence.getDataXYZAsByte(0, 0);
    getSeedHistograms();
    float val, maxVal = 0;
    if (createProbas)
      terminalSequence.beginUpdate();
    for (z = 0; z < sz; z++) {
      float[] terminalData = null;
      if (createProbas) {
        terminalSequence.setImage(0, z,
            new IcyBufferedImage(sx, sy, 1, DataType.FLOAT));
        terminalData = terminalSequence.getDataXYAsFloat(0, z, 0);
      }
      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          currPos = z * sxy + y * sx + x;

          probaSrc = new Double(
              -Math.log(getProba(0, TypeUtil.unsign(grayData[z][y * sx + x]))))
                  .floatValue();
          // System.out.print("proba (" + probaSrc);
          probaSnk = new Double(
              -Math.log(getProba(1, TypeUtil.unsign(grayData[z][y * sx + x]))))
                  .floatValue();
          // System.out.println(", " + probaSnk+")");
          val = graphCut.setTerminalWeights(currPos, probaSrc, probaSnk, lambda,
              false);
          maxVal = Math.max(maxVal, Math.abs(val));
          if (createProbas)
            terminalData[y * sx + x] = val;
        }
      }
    }

    maxVal += 1;
    float[][] terminalData = null;
    if (createProbas)
      terminalData = terminalSequence.getDataXYZAsFloat(0, 0);
    SequenceDataIterator sDI;
    for (i = 0; i < 2; i++) {
      sDI = new SequenceDataIterator(treatedSequence, seeds.get(i));
      while (!sDI.done()) {
        currPos = sDI.getPositionZ() * sxy + sDI.getPositionY() * sx
            + sDI.getPositionX();
        val = graphCut.setTerminalWeights(currPos, (float) (i), (float) (1 - i),
            maxVal, false);
        if (createProbas)
          terminalData[sDI.getPositionZ()][sDI.getPositionY() * sx
              + sDI.getPositionX()] = val;
        sDI.next();
      }
    }
    if (createProbas) {
      terminalSequence.dataChanged();
      terminalSequence.endUpdate();
    }

  }

  private double getProba(int seed, int level) {
    double count = 0;
    // double div = 0;
    count = seedHistograms[seed].getBin(level).getCount();
    /*
     * for (int i = -edgeStdDev; i <= edgeStdDev; i++) { if (level + i >= 0 &&
     * level + i < seedHistograms[seed].getNbBins()) { count +=
     * seedHistograms[seed].getBin(level+i).getCount(); div += 1; } } count =
     * count/div;
     */
    return count / seedSizes[seed];
  }

  /**
   * Computes the gray level variance of the treated sequence
   * 
   * @return gray level variance
   */
  private double computeVariance() {
    int z, xy;

    byte[][] seqData = graySequence.getDataXYZAsByte(0, 0);
    int sxy = sx * sy;
    double mean = 0;
    for (z = 0; z < sz; z++) {
      for (xy = 0; xy < sxy; xy++) {
        mean += seqData[z][xy];
      }
    }
    mean /= sxy * sz;
    System.out.println(mean);

    double variance = 0;
    double val;
    for (z = 0; z < sz; z++) {
      for (xy = 0; xy < sxy; xy++) {
        val = mean - seqData[z][xy];
        variance += val * val;
      }
    }
    variance /= sxy * sz;
    System.out.println(variance);
    return variance;
  }

  private double getEdgeLikelihood(int[] values1, int[] values2, int dx, int dy,
      int dz) {
    double diff = 0;
    for (int c = 0; c < sc; c++) {
      diff = Math.max(diff, Math.abs(values2[c] - values1[c]));
    }
    double val = Math.exp(-(diff * diff) / (2.0 * edgeVariance))
        / Math.sqrt(dx * dx + dy * dy + dz * dz);
    return val;
  }

  private Histogram[] getSeedHistograms() throws BadHistogramParameters {
    int i;

    seedSizes = new double[seeds.size()];
    for (i = 0; i < seedSizes.length; i++) {
      seedSizes[i] = seeds.get(i).getNumberOfPoints();
    }

    seedHistograms = new Histogram[seeds.size()];
    for (i = 0; i < seeds.size(); i++) {
      seedHistograms[i] =
          Histogram.compute(graySequence, seeds.get(i), true, 256, 0, 255);
    }
    return seedHistograms;
  }

  @Override
  public void executeSegmentation() {
    int x, y, z, r;
    int sxy = sx * sy;

    double maxFlow = graphCut.computeMaximumFlow(false, null);
    System.out.printf("Max flow = %f%n", maxFlow);

    List<ROI3DArea> resROI = new ArrayList<ROI3DArea>();
    for (r = 0; r < 2; r++) {
      resROI.add(new ROI3DArea());
    }

    for (z = 0; z < sz; z++) {

      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          resROI.get((graphCut
              .getTerminal(z * sxy + y * sx + x) == Terminal.BACKGROUND)? 0: 1)
              .addPoint(x, y, z);
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
      try {
        for (int z = 0; z < sz; z++)
          resultSegmentation.setImage(0, z,
              new IcyBufferedImage(sx, sy, 1, DataType.UBYTE));

        // set value from ROI(s)
        byte val = (int)0;
        for (ROI roi: resultROIs) {
          if (!roi.getBounds5D().isEmpty())
            DataIteratorUtil
                .set(new SequenceDataIterator(resultSegmentation, roi), val);
          val = (byte)DataType.UBYTE_MAX_VALUE;
        }

        // notify data changed
        resultSegmentation.dataChanged();
      }
      finally {
        resultSegmentation.endUpdate();
      }
      resultSegmentation
          .setName(inSequence.getName() + "(var" + String.format("%.2f", edgeVariance) + ", lambda"
              + lambda + ", " + ((use8Connected)? "8": "4") + "-conn)");
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
    result.setName(inSequence.getName() + "Segmentation(var" + String.format("%.2f", edgeVariance)
        + ", lambda" + lambda + ", " + ((use8Connected)? "8": "4") + "-conn)");
    return result;
  }

  public Sequence getTreatedSequence() {
    return this.treatedSequence;
  }

  public Sequence getGradientSequence() {
    return gradientSequence;
  }

  public Sequence getTerminalProbabilitiesSequence() {
    return terminalSequence;
  }

}
