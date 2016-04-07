/**
 * 
 */
package algorithms.danyfel80.segmentation.graphcut.old;

import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import icy.type.point.Point5D;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;

import javax.vecmath.Point3i;

import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

import com.google.common.collect.Lists;

import algorithms.danyfel80.segmentation.SegmentationAlgorithm;


/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class GraphCutSegmentationHashMap extends SegmentationAlgorithm {

  private Sequence inSequence;
  private Sequence grayInSequence;
  private Sequence segmentation;
  private List<ROI> seeds;

  private int sizeX;
  private int sizeY;
  private int sizeZ;

  private final double lambda;
  private final double K;
  private DefaultDirectedWeightedGraph<Point3i, DefaultWeightedEdge> graph;

  private Color[] colors;

  private static final Point3i SOURCE = new Point3i(-1, -1, -1);
  private static final Point3i TARGET = new Point3i(-2, -2, -2);
  private Map<Point3i, Integer> trees;
  private Map<Point3i, Point3i> parents;
  private Queue<Point3i> activeNodes;
  private Queue<Point3i> orphanNodes;

  public GraphCutSegmentationHashMap(Sequence inSequence, double lambda, List<ROI> seeds) {
    this.K = 27.0;
    this.lambda = lambda;
    this.inSequence = inSequence;
    this.seeds = seeds;
    prepareGraph(seeds);
  }

  /* (non-Javadoc)
   * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#prepareGraph()
   */
  @Override
  protected void prepareGraph(List<ROI> seeds) {
    grayInSequence = SequenceUtil.toGray(inSequence);
    grayInSequence = SequenceUtil.convertToType(grayInSequence, DataType.DOUBLE, false);
    double[][][] grayData = grayInSequence.getDataXYCZAsDouble(0);
    segmentation = new Sequence();

    sizeZ = grayInSequence.getSizeZ();
    sizeY = grayInSequence.getSizeY();
    sizeX = grayInSequence.getSizeX();

    graph = new DefaultDirectedWeightedGraph<>(DefaultWeightedEdge.class);

    // Add vertices to graph, compute mean intensity, create segmentation object
    double meanI = 0.0;
    segmentation.beginUpdate();
    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          graph.addVertex(new Point3i(x, y, z)); // add node (x,y,z)
          meanI += grayData[z][0][x + y*sizeX];
        }
      }
      segmentation.setImage(0, z, new IcyBufferedImage(sizeX, sizeY,3, DataType.DOUBLE));
    }
    segmentation.endUpdate();
    graph.addVertex(SOURCE); // add source supplementary node
    graph.addVertex(TARGET); // add target supplementary node

    meanI /= (double)(sizeX*sizeY*sizeZ);

    // get variance
    double varianceI = 0.0;
    double diff = 0;
    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          diff = grayData[z][0][x + y*sizeX] - meanI;
          varianceI += diff*diff;
        }
      }
    }
    varianceI /= (double)(sizeX*sizeY*sizeZ);

    // Add n-links
    for (int z = 0; z < sizeZ; z++) {
      for (int y = 0; y < sizeY; y++) {
        for (int x = 0; x < sizeX; x++) {

          for (int dz = -1; dz <= 1; dz++) {
            for (int dy = -1; dy <= 1; dy++) {
              for (int dx = -1; dx <= 1; dx++) {

                if (dz != 0 || dy != 0 || dx != 0) {
                  if (z+dz >= 0 && z+dz < sizeZ &&
                      y+dy >= 0 && y+dy < sizeY &&
                      x+dx >= 0 && x+dx < sizeX) {
                    Point3i p = new Point3i(x, y, z);
                    Point3i q = new Point3i(x+dx, y+dy, z+dz);
                    DefaultWeightedEdge e1 = graph.addEdge(p, q);
                    DefaultWeightedEdge e2 = graph.addEdge(q, p);
                    if (e1 != null && e2 != null){
                      diff = grayData[z][0][x+y*sizeX] - grayData[(z+dz)][0][(x+dx)+(y+dy)*sizeX];
                      Double weight = Math.exp(-(diff*diff)/(2.0*varianceI)) * (1.0/Math.sqrt(dx*dx+dy*dy+dz*dz));
                      graph.setEdgeWeight(e1, weight);
                      graph.setEdgeWeight(e2, weight);
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

  /* (non-Javadoc)
   * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#executeSegmentation(java.util.List)
   */
  @Override
  public void executeSegmentation() {
    // compute t-links
    List<Color> labels = new ArrayList<>();
    for (ROI roi : seeds) {
      labels.add(roi.getColor());
    }
    colors = labels.toArray(new Color[2]);
    paintSeeds(seeds);
    double[][] histograms = getClassesHistograms();
    addTLinks(histograms);

    // max-flow algorithm

    // 1. init
    trees = new HashMap<Point3i, Integer>();
    parents = new HashMap<Point3i, Point3i>();
    activeNodes = new LinkedList<Point3i>();
    orphanNodes = new LinkedList<Point3i>();


    trees.put(SOURCE, 0);
    trees.put(TARGET, 1);
    activeNodes.add(SOURCE);
    activeNodes.add(TARGET);

    while (true) {
      //System.out.println("SizeTrees = " + trees.size());
      // 2. Grow
      List<Point3i> augPath = growTrees();
      if (augPath.isEmpty())
        break;
      //System.out.println("Path: " + augPath);
      // 3. Augment
      augmentPath(augPath);
      //System.out.println("Augmented, orphans = " + orphanNodes);
      // 4. Adopt
      adoptOrphanNodes();
      //System.out.println("Adopted");
      //System.out.flush();
    }

    paintSegmentation();
  }

  private void paintSeeds(List<ROI> seeds) {
    segmentation.beginUpdate();
    double[][][] segData = segmentation.getDataXYCZAsDouble(0);
    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          if (seeds.get(0).contains(new Point5D.Double(x,y,z,0,0)))
            segData[z][0][x*sizeY + y] = 1;
          else if (seeds.get(1).contains(new Point5D.Double(x,y,z,0,0)))
            segData[z][0][x*sizeY + y] = 2;
        }
      }
    }
    segmentation.dataChanged();
    segmentation.endUpdate();
  }

  /**
   * Computes the histogram for each of the seed ROIs.
   * @return Histogram array (one histogram for each ROI).
   */
  private double[][] getClassesHistograms() {
    int [] amount = new int[2];
    double[][] histograms = new double[2][256];

    double [][][] segData = segmentation.getDataXYCZAsDouble(0); 
    double [][][] imageData = grayInSequence.getDataXYCZAsDouble(0);

    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          double v = segData[z][0][x + y*sizeX];
          if (v > 0) {
            int p = ((int)imageData[z][0][x + y*sizeX]);
            histograms[(int)v-1][p] += 1;
            amount[(int)v-1]++;
          }
        }
      }
    }
    for (int i = 0; i < 256; i++) {
      histograms[0][i] /= (double)amount[0];
      histograms[1][i] /= (double)amount[1];
      //System.out.printf("[%d]=%f\n", i, histograms[0][i]);
    }

    return histograms;
  }

  private void addTLinks(double[][] histograms) {
    double [][][] grayData = grayInSequence.getDataXYCZAsDouble(0);
    double [][][] segData = segmentation.getDataXYCZAsDouble(0);

    for (int z = 0; z < sizeZ; z++) {
      for (int x = 0; x < sizeX; x++) {
        for (int y = 0; y < sizeY; y++) {
          Point3i pos = new Point3i(x, y, z);
          DefaultWeightedEdge eS1 = graph.addEdge(pos, SOURCE);
          DefaultWeightedEdge eS2 = graph.addEdge(SOURCE, pos);
          DefaultWeightedEdge eT1 = graph.addEdge(pos, TARGET);
          DefaultWeightedEdge eT2 = graph.addEdge(TARGET, pos);
          Double weightS = 0.0;
          Double weightT = 0.0;
          if (segData[z][0][x*sizeY + y] > 0.0) { // if seed
            weightS = (segData[z][0][x*sizeY + y] == 2.0)? K: 0.0;
            weightT = (segData[z][0][x*sizeY + y] == 1.0)? K: 0.0;

          } else { // if not seed, then use histograms
            int p = (int)grayData[z][0][x + y*sizeX];
            weightS = lambda * (-Math.log(histograms[0][p]));
            weightT = lambda * (-Math.log(histograms[1][p]));
            if (weightS > K) weightS = K;
            if (weightT > K) weightT = K;
          }
          graph.setEdgeWeight(eS1, weightS);
          graph.setEdgeWeight(eS2, weightS);
          graph.setEdgeWeight(eT1, weightT);
          graph.setEdgeWeight(eT2, weightT);
        }
      }
    }
  }

  private List<Point3i> growTrees() {
    while (!activeNodes.isEmpty()) {
      //System.out.println("Active nodes: " + activeNodes.size());
      Point3i p = activeNodes.element();
      Set<DefaultWeightedEdge> neighbors = graph.outgoingEdgesOf(p);
      for (DefaultWeightedEdge pq : neighbors) { // search all neighbors with capacity > 0
        Point3i q = graph.getEdgeTarget(pq);
        if (getEdgeCapacity(p, q) > 0.0) {
          Integer treeQ = trees.get(q);
          Integer treeP = trees.get(p);
          if (treeQ == null) {
            trees.put(q, treeP);
            parents.put(q, p);
            activeNodes.add(q);
          } 
          if (treeQ != null && !treeQ.equals(treeP)) { // Found Augmenting Path!!
            List<Point3i> path = new ArrayList<Point3i>();
            Point3i node = p;
            while (node != null) {
              path.add(node);
              node = parents.get(node);
            }
            path = Lists.reverse(path);
            node = q;
            while (node != null) {
              path.add(node);
              node = parents.get(node);
            }
            if (path.get(0).equals(TARGET)) {
              path = Lists.reverse(path);
            }
            return path;
          }
        }
      }
      activeNodes.remove();
    }
    return new ArrayList<Point3i>();
  }

  private void augmentPath(List<Point3i> augPath) {
    List<DefaultWeightedEdge> pathEdges = new ArrayList<DefaultWeightedEdge>(augPath.size()-1);
    List<DefaultWeightedEdge> revPathEdges = new ArrayList<DefaultWeightedEdge>(augPath.size()-1);

    // Find bottleneck capacity
    double maxCapacity = Double.MAX_VALUE;
    for (int i = 1; i < augPath.size(); i++) {
      DefaultWeightedEdge e = graph.getEdge(augPath.get(i-1), augPath.get(i));
      pathEdges.add(e);
      revPathEdges.add(graph.getEdge(augPath.get(i), augPath.get(i-1)));
      maxCapacity = Math.min(maxCapacity, graph.getEdgeWeight(e));
    }

    // Update residual graph by pushing maxCapacity throw augPath.
    List<Double> newWeights = new ArrayList<>(pathEdges.size());
    for (int i = 0; i < pathEdges.size(); i++) {
      DefaultWeightedEdge ei = pathEdges.get(i);
      DefaultWeightedEdge eRevi = revPathEdges.get(i);
      double newWeight = graph.getEdgeWeight(ei) - maxCapacity;
      double newRevWeight = graph.getEdgeWeight(eRevi) + maxCapacity;
      graph.setEdgeWeight(ei, newWeight);
      graph.setEdgeWeight(eRevi, newRevWeight);
      newWeights.add(newWeight);
    }

    // Find saturated edges and create orphan nodes
    for (int i = 0; i < pathEdges.size(); i++) {
      if (newWeights.get(i) == 0.0) {
        Point3i p = augPath.get(i);
        Point3i q = augPath.get(i+1);
        Integer treeP = trees.get(p);
        Integer treeQ = trees.get(q);
        if (treeP.equals(0) && treeP.equals(treeQ)) {
          parents.remove(q);
          orphanNodes.add(q);
        }
        if (treeP.equals(1) && treeP.equals(treeQ)) {
          parents.remove(p);
          orphanNodes.add(p);
        }
      }
    }
  }

  private void adoptOrphanNodes() {
    //boolean debug = true;
    //Point3i pDeb = new Point3i(23, 8, 0);
    while (!orphanNodes.isEmpty()) {
      Point3i p = orphanNodes.remove();
      //if (p.equals(pDeb)) debug = true;
      //if (debug) System.out.println(orphanNodes.size());
      Set<DefaultWeightedEdge> neighborEdges = graph.outgoingEdgesOf(p);

      boolean foundParent = false;
      for (DefaultWeightedEdge neighborEdge : neighborEdges) {
        Point3i q = graph.getEdgeTarget(neighborEdge);
        Integer treeP = trees.get(p);
        Integer treeQ = trees.get(q);
        if ((treeP != null && treeQ != null && treeQ.equals(treeP)) && // 
            getEdgeCapacity(q, p) > 0.0) {
          Point3i origin = q;
          Point3i pOrigin = parents.get(origin);
          while(pOrigin != null && !pOrigin.equals(p)){
            //System.out.println("Parent(" + origin + ")=" + pOrigin);
            origin = pOrigin;
            pOrigin = parents.get(origin);
          }
          if (origin.equals(SOURCE) || origin.equals(TARGET)) {
            parents.put(p, q);
            foundParent = true;
          }
        }
      }

      if (!foundParent) {
        for (DefaultWeightedEdge neighborEdge : neighborEdges) {
          Point3i q = graph.getEdgeTarget(neighborEdge);
          Integer treeP = trees.get(p);
          Integer treeQ = trees.get(q);
          if (treeP != null && treeQ != null && treeQ.equals(treeP)) {
            if (getEdgeCapacity(q, p) > 0.0)
              activeNodes.add(q);
            Point3i parQ = parents.get(q);
            if (parQ != null && parQ.equals(p)) {
              orphanNodes.add(q);
              parents.remove(q);
            }
          }
        }
        trees.remove(p);
        while(activeNodes.remove(p));
      }
    }
  }

  private Double getEdgeCapacity(Point3i p, Point3i q) {
    Integer treeVal = trees.get(p);
    if (treeVal == null) // if p is free
      return null;
    else if (treeVal == 0) // if p belongs to S
      return graph.getEdgeWeight(graph.getEdge(p, q));
    else // if p belongs to T
      return graph.getEdgeWeight(graph.getEdge(q, p));
  }

  private void paintSegmentation() {
    segmentation.beginUpdate();
    double[][][] segData = segmentation.getDataXYCZAsDouble(0);
    Set<Entry<Point3i, Integer>> treeEntries = trees.entrySet();
    for (Entry<Point3i, Integer> entry : treeEntries) {
      Point3i coord = entry.getKey();
      if (!coord.equals(SOURCE) && !coord.equals(TARGET)) {
        segData[coord.z][0][coord.x + coord.y*sizeX] = colors[(int)entry.getValue()].getRed();
        segData[coord.z][1][coord.x + coord.y*sizeX] = colors[(int)entry.getValue()].getGreen();
        segData[coord.z][2][coord.x + coord.y*sizeX] = colors[(int)entry.getValue()].getBlue();
      }
    }
    segmentation.dataChanged();
    segmentation.endUpdate();
  }

  @Override
  public Sequence getSegmentationSequence() {
    // TODO Auto-generated method stub
    System.err.println("Method not yet implemented.");
    return null;
  }

  @Override
  public Collection<? extends ROI> getSegmentationROIs() {
    // TODO Auto-generated method stub
    System.err.println("Method not yet implemented.");
    return null;
  }

  @Override
  public Sequence getSegmentationSequenceWithROIs() {
    return segmentation;
  }

}
