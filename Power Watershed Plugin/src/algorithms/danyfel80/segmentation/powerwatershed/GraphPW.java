package algorithms.danyfel80.segmentation.powerwatershed;

import java.awt.Color;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Deque;
import java.util.List;

import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceDataIterator;
import icy.type.DataType;
import icy.type.TypeUtil;
import plugins.kernel.roi.roi3d.ROI3DArea;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class GraphPW {
  private int N; // #nodes
  private int M; // #edges
  private Sequence seq; // node data
  private int sx, sy, sz; // sequence size data
  private List<ROI> seeds; // seed data

  private boolean useGeo;
  private boolean useQuickSort;
  private int[][] edges; // edges source-target
  int numSeeds; // #seeds
  private int[] indexSeeds; // seed to node index
  private int[] indexLabels; // seed to label index
  int numLabels; // #labels
  Color[] colorLabels; // color to label index
  String[] nameLabels; // name of each label

  private int[] weights; // edge weights
  private int[] normalWeights; // edge normalized weights

  public static final double EPSILON = 1e-6;
  public static final int MAX_SIZE_PLATEAU = 1000000;

  private Sequence seqResult;

  /**
   * Constructor for the power watershed segmentation algorithm
   * @param N Amount of nodes
   * @param M Amount of edges
   * @param seq Sequence where the segmentation is performed
   * @param seeds Seed ROIS. One ROI per label.
   */
  public GraphPW(Sequence seq, List<ROI> seeds, boolean useGeo, boolean useQuickSort) { 

    this.seq = seq;
    this.seeds = seeds;
    this.useGeo = useGeo;

    sx = seq.getSizeX(); // Size of sequence in x
    sy = seq.getSizeY(); // Size of sequence in y
    sz = seq.getSizeZ(); // Size of sequence in z

    N = sx*sy*sz; // # nodes
    M = 3*sx*sy*sz - sx*sy - sx*sz - sy*sz; // # edges

    edges = new int[2][]; // Edges (pairs of source->target)
    for (int i = 0; i < 2; i++) {
      edges[i] = new int[M];
    }
    indexSeeds = new int[N]; // Seed index (position in image of each seed node)
    indexLabels = new int[N]; // Label index (label of each node)

    seqResult = null; 
  }

  /**
   * Creates the index of seeds and captures the seeds' colors.
   */
  public void SetupSeeds() {
    numLabels = seeds.size();
    colorLabels = new Color[numLabels];
    nameLabels = new String[numLabels];

    int x, y, z;
    numSeeds = 0;
    SequenceDataIterator iter;
    for (int i = 0; i < numLabels; i++) {
      colorLabels[i] = seeds.get(i).getColor();
      nameLabels[i] = seeds.get(i).getName();
      System.out.printf("color %d for " + nameLabels[i] + " = " + colorLabels[i] + "\n", i);
      iter = new SequenceDataIterator(seq, seeds.get(i));
      iter.reset();
      while (!iter.done()) {
        x = iter.getPositionX();
        y = iter.getPositionY();
        z = iter.getPositionZ();

        indexSeeds[numSeeds] = z*sx*sy + y*sx + x;
        indexLabels[numSeeds] = i;
        numSeeds++;
        iter.next();
      }
    }
  }

  /**
   * Sets the connections between nodes.
   */
  public void SetupEdges() {
    int sxy = sx*sy; // size of level
    int x, y, z, l, M = 0; 
    for (z = 0; z < sz; z++) {
      // (x,y,z) -> (x,y+1,z)
      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          if (y < sy-1) {
            edges[0][M] = x + y*sx + z*sxy;
            edges[1][M] = x + (y+1)*sx + z*sxy;
            M++;
          }
        }
      }
      // (x,y,z) -> (x+1,y,z)
      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          if (x < sx-1) {
            edges[0][M] = x + y*sx + z*sxy;
            edges[1][M] = (x+1) + y*sx + z*sxy;
            M++;
          }
        }
      }
      // (x,y,z) -> (x,y,z+1)
      if (z < sz-1) {
        for (l = z*sxy; l < (z+1)*sxy; l++) { // every node in level z
          edges[0][M] = l;
          edges[1][M] = l + sxy;
          M++;
        }
      }
    }
  }

  /**
   * Adds the weights of each edge and performs geodesic reconstruction.
   */
  public void calculateWeights() {
    byte[][][] seqData = seq.getDataXYCZAsByte(0);
    int i, j, xy1, z1, xy2, z2, sxy = sx*sy;
    
    int[] vals = new int[seq.getSizeC()];
    normalWeights = new int[M];
    for (i = 0; i < M; i++) {
      z1 = edges[0][i]/sxy;
      xy1 = edges[0][i]%sxy;
      z2 = edges[1][i]/sxy;
      xy2 = edges[1][i]%sxy;
      for (j = 0; j < seq.getSizeC(); j++) {
        vals[j] = Math.abs(TypeUtil.unsign(seqData[z1][j][xy1]) - TypeUtil.unsign(seqData[z2][j][xy2]));
        if (j == 0) {
          normalWeights[i] = 255 - vals[j];
        } else {
          normalWeights[i] = Math.min(normalWeights[i], 255 - vals[j]);
        }
      }
    }

    
    int numNeighbors = (sz == 1)? 4: 6;
    int n;

    int[] seedsFunction = new int[M];
    for (i = 0; i < numSeeds; i++) {
      for (j = 1; j <= numNeighbors; j++) {
        n = GraphUtils.getNeighborNodeEdge(indexSeeds[i], j, sx, sy, sz);
        if (n != -1) {
          seedsFunction[n] = normalWeights[n];
        }
      }
    }
    
    weights = new int[M];
    // if not dilated the result of segmentation is voronoi because weights are all 0.
    DilateReconstruction.reconstruct(seedsFunction, normalWeights, weights, edges, sx, sy, sz, M, useQuickSort);
  }
  
  public Sequence getSequenceGradient() {
    int k, l, x, y, z, n;
    double w;
    int numNeighbors = (sz == 1)? 4: 6;
    Sequence result = new Sequence(seq.getName() + "Gradient");
    result.beginUpdate();
    for (z = 0; z < sz; z++) {
      result.setImage(0, z, new IcyBufferedImage(sx, sy, 1, DataType.DOUBLE));
    }
    
    double[][] resData = result.getDataXYZAsDouble(0, 0);
    for (z = 0; z < sz; z++) {
      for (y = 0; y < sy; y++) {
        for (x = 0; x < sx; x++) {
          l = 0;
          w = 0;
          for (k = 1; k <= numNeighbors; k++) {
            n = GraphUtils.getNeighborNodeEdge(x + y*sx + z*sy*sx, k, sx, sy, sz);
            if (n != -1) {
              w += weights[n];
              l++;
            }
          }
          resData[z][x + y*sx] = ((2.0*w)/(double)l);
        }
      }
    }
    result.dataChanged();
    result.endUpdate();
    
    return result;
  }

  /**
   * Executes the power watershed segmentation method.
   * @param useGeo Make geodesic reconstruction.
   */
  public void executeWatershed() {
    if(useGeo) {
      seqResult = executePWQ2(edges, weights, weights, 255, indexSeeds, indexLabels, numSeeds, sx, sy, sz, numLabels);
    }
    else {
      seqResult = executePWQ2(edges, weights, normalWeights, 255, indexSeeds, indexLabels, numSeeds, sx, sy, sz, numLabels);
    }
  }

  /**
   * Computes the result sequence of the energy minimization : min_x lim_p->inf sum_e_of_E w_{ij}^p |x_i-x_j|^2 
   * @param edges Array of node indexes composing edges
   * @param weights Reconstructed weights
   * @param normalWeights Original weights
   * @param maxWeight Maximum weight value
   * @param indexSeeds Array of seeded nodes' indexes
   * @param indexLabels Label values on the seeded nodes
   * @param numSeeds Amount of seeded nodes
   * @param sx X-coordinate size of the sequence
   * @param sy Y-coordinate size of the sequence
   * @param sz Z-coordinate size of the sequence
   * @param numLabels Number of different labels
   * @return The result x of the energy function minimization
   */
  private Sequence executePWQ2(int[][] edges, int[] weights,
      int[] normalWeights, int maxWeight, int[] indexSeeds, int[] indexLabels,
      int numSeeds, int sx, int sy, int sz, int numLabels) {

    int[] fathers = new int[N];
    int[] localSeeds = new int[N];
    int[] sortedEdges = new int[M];
    int[] sortedNormalEdges = new int[M];
    boolean[] edgeIndicators = new boolean[M];
    boolean[] plateauIndicators = new boolean[M];
    int[] plateauValueIndicators = new int[N];
    int[] plateauVerticesList = new int[N];
    int[] ranks = new int[N];
    Deque<Integer> plateauStack = new ArrayDeque<Integer>(M);  
    Deque<Integer> LCPStack = new ArrayDeque<Integer>(M);

    int i, j, k, x, y, e1, e2, re1, re2, p, xr;
    int argMax;
    int numVertices, numEdges, numNormalEdges, edgeMax, normalEdgeMax, maxW;
    int numNeighborEdges = (sz > 1)? 12: 6;

    boolean success = false, differentSeeds;

    double val;

    double[][] proba = new double[numLabels - 1][];
    for (i = 0; i < numLabels - 1; i++) {
      proba[i] = new double[N];
      for (j = 0; j < N; j++) {
        proba[i][j] = -1;
      }
    }

    int[][] edgesLCP = new int[2][];
    for (i = 0; i < 2; i++) {
      edgesLCP[i] = new int[M];
    }

    for (i = 0; i < numSeeds; i++) {
      for (j = 0; j < numLabels - 1; j++) {
        if (indexLabels[i] == j+1) {
          proba[j][indexSeeds[i]] = 1;
        } else {
          proba[j][indexSeeds[i]] = 0;
        }
      }
    }

    for ( i = 0; i < N; i++) {
      fathers[i] = i;
    }

    double[][] localLabels = new double[numLabels - 1][];
    for (i = 0; i < numLabels - 1; i++) {
      localLabels[i] = new double[N];
    }

    int[] sortedWeights = new int[M];
    for (i = 0; i < M; i++) {
      sortedWeights[i] = weights[i];
      sortedEdges[i] = i;
    }

    if (useQuickSort) {
      GraphUtils.binSortDec(sortedWeights, sortedEdges, M, maxWeight + 1);
    } else {
      GraphUtils.quickStochasticSortDec(sortedWeights, sortedEdges, 0, M - 1);
    }

    int edgeCount = 0;
    int normalEdgeCount = 0;

    // Main loop
    while (edgeCount < M) {
      do {
        edgeMax = sortedEdges[edgeCount++];
        if (edgeCount >= M) {
          break;
        }
      } while(edgeIndicators[edgeMax]);

      if (edgeCount >= M) {
        break;
      }

      // 1. Computing the edges of the plateau LCP linked to the edge edgeMax
      plateauStack.addFirst(edgeMax);
      plateauIndicators[edgeMax] = true;
      edgeIndicators[edgeMax] = true;
      LCPStack.addFirst(edgeMax);
      numVertices = 0;
      numEdges = 0;
      maxW = weights[edgeMax];

      // 2. Putting the edges and vertices of the plateau into arrays
      while(!plateauStack.isEmpty()) {
        x = plateauStack.removeFirst();
        e1 = edges[0][x];
        e2 = edges[1][x];
        re1 = GraphUtils.findElement(e1, fathers);
        re2 = GraphUtils.findElement(e2, fathers);

        if (proba[0][re1] < 0 || proba[0][re2] < 0) {
          if (plateauValueIndicators[e1] == 0) {
            plateauVerticesList[numVertices++] = e1;
            plateauValueIndicators[e1] = 1;
          }

          if (plateauValueIndicators[e2] == 0) {
            plateauVerticesList[numVertices++] = e2;
            plateauValueIndicators[e2] = 1;
          }

          edgesLCP[0][numEdges] = e1;
          edgesLCP[1][numEdges] = e2;
          sortedNormalEdges[numEdges++] = x;
        }

        for (k = 1; k <= numNeighborEdges; k++) {

          y = ((sz > 1)? 
              GraphUtils.getNeighborEdge3D(e1, e2, x, k, sx, sy, sz): 
                GraphUtils.getNeighborEdge(x, k, sx, sy, sz));

          if (y != -1) {
            if ((!plateauIndicators[y]) && (weights[y] == maxW)) {
              plateauIndicators[y] = true;
              plateauStack.addFirst(y);
              LCPStack.addFirst(y);
              edgeIndicators[y] = true;
            }
          }
        }
      }

      for (j = 0; j < numVertices; j++) {
        plateauValueIndicators[plateauVerticesList[j]] = 0;
      }

      for (int l : LCPStack) {
        plateauIndicators[l] = false;
      }

      // 3. If maxEdge belongs to a plateau
      if (numEdges > 0) {
        // 4. Evaluate if there are differents seeds on the plateau
        p = 0;
        differentSeeds = false;

        for (i = 0; i < numLabels - 1; i++) {
          val = -0.5;
          for (j = 0; j < numVertices; j++) {
            x = plateauVerticesList[j];
            xr = GraphUtils.findElement(x, fathers);
            if (Math.abs(proba[i][xr] - val) > EPSILON && proba[i][xr] >= 0) {
              p++;
              val = proba[i][xr];
            }
          }

          if (p >= 2) {
            differentSeeds = true;
            break;
          } else {
            p = 0;
          }
        }

        if (differentSeeds) {
          // 5. Sort the edges of the plateau according to their normal weight
          for (k = 0; k < numEdges; k++) {
            sortedWeights[k] = normalWeights[sortedNormalEdges[k]];
          }
          
          
          if (useQuickSort) {
            GraphUtils.binSortDec(sortedWeights, sortedNormalEdges, numEdges, maxWeight + 1);
          } else {
            GraphUtils.quickStochasticSortDec(sortedEdges, sortedNormalEdges, 0, numEdges - 1);
          }
          
          // Merge nodes for edges of real max weight
          numVertices = 0;
          numNormalEdges = 0;
          for (normalEdgeCount = 0; normalEdgeCount < numEdges; normalEdgeCount++) {
            normalEdgeMax = sortedNormalEdges[normalEdgeCount];
            e1 = edges[0][normalEdgeMax];
            e2 = edges[1][normalEdgeMax];
            if (normalWeights[normalEdgeMax] != maxW) {
              mergeNode(e1, e2, ranks, fathers, proba, numLabels);
            } else {
              re1 = GraphUtils.findElement(e1, fathers);
              re2 = GraphUtils.findElement(e2, fathers);
              if (re1 != re2 && (proba[0][re1] < 0 || proba[0][re2] < 0)) {
                if (plateauValueIndicators[re1] == 0) {
                  plateauVerticesList[numVertices++] = re1;
                  plateauValueIndicators[re1] = 1;
                }

                if (plateauValueIndicators[re2] == 0) {
                  plateauVerticesList[numVertices++] = re2;
                  plateauValueIndicators[re2] = 1;
                }

                edgesLCP[0][numNormalEdges] = re1;
                edgesLCP[1][numNormalEdges] = re2;
                numNormalEdges++;
              }
            }
          }

          for (i = 0; i < numLabels - 1; i++) {
            k = 0;
            for (j = 0; j < numVertices; j++) {
              xr = plateauVerticesList[j];
              if (proba[i][xr] >= 0) {
                localLabels[i][k] = proba[i][xr];
                localSeeds[k++] = xr;
              }
            }
          }

          // 6. Execute Random Walker on plateaus
          if (numVertices < MAX_SIZE_PLATEAU) {
            success = RandomWalker.ExecuteRandomWalker(edgesLCP, numNormalEdges, plateauVerticesList, plateauValueIndicators, numVertices, localSeeds, localLabels, k, numLabels, proba);
          }

          if (numVertices >= MAX_SIZE_PLATEAU || !success) {
            System.out.printf("Plateau of a big size (%d vertices, %d edges) the RW is not performed on it.\n", numVertices, numNormalEdges);
            for (j = 0; j < numNormalEdges; j++) {
              e1 = edgesLCP[0][j];
              e2 = edgesLCP[1][j];
              mergeNode(e1, e2, ranks, fathers, proba, numLabels);
            }
          }

          for (j = 0; j < numVertices; j++) {
            plateauValueIndicators[plateauVerticesList[j]] = 0;
          }
        } else { // if different seeds = false
          // 7. Merge nodes for edges of max weight
          for (j = 0; j < numEdges; j++) {
            e1 = edgesLCP[0][j];
            e2 = edgesLCP[1][j];
            mergeNode(e1, e2, ranks, fathers, proba, numLabels);
          }
        }
      }

      LCPStack.clear();
    } // End main loop

    // Building the final proba map (find the root vertex of each tree)
    for (i = 0; i < N; i++) {
      j = i;
      xr = i;

      while (fathers[i] != i) {
        i = xr;
        xr = fathers[i];
      }

      for (k = 0; k < numLabels - 1; k++) {
        proba[k][j] = proba[k][i];
      }

      i = j;
    }

    // Write results
    seqResult = new Sequence(this.seq.getName());
    seqResult.beginUpdate();
    for (int z = 0; z < sz; z++) {
      seqResult.setImage(0, z, new IcyBufferedImage(sx, sy, 1, DataType.INT));
    }


    int[][] resultData = seqResult.getDataXYZAsInt(0, 0);
    double maxi;
    int sxy = sx*sy;

    for (j = 0; j < N; j++) {
      maxi = 0;
      argMax = 0;
      val = 1;
      for (k = 0; k < numLabels - 1; k++) {
        if (proba[k][j] > maxi) {
          maxi = proba[k][j];
          argMax = k;
        }
        val = val - proba[k][j];
      }

      if (val > maxi) {
        argMax = k;
      }

      resultData[j/sxy][j%sxy] = argMax;//(argMax*255)/(numLabels - 1);
    }
    seqResult.dataChanged();
    seqResult.endUpdate();
    return seqResult;
  }

  /**
   * Updates the result(proba), ranks and fathers arrays when 2 nodes are merged
   * @param e1 Index of node 1
   * @param e2 Index of node 2
   * @param ranks Array needed for union-find efficiency
   * @param fathers Array for storing roots of merged nodes trees
   * @param proba Array for storing the result x
   * @param numLabels Amount of labels
   */
  private void mergeNode(int e1, int e2, int[] ranks, int[] fathers,
      double[][] proba, int numLabels) {
    int k, re1, re2;
    re1 = GraphUtils.findElement(e1, fathers);
    re2 = GraphUtils.findElement(e2, fathers);

    if(re1 != re2 && !(proba[0][re1] >= 0 && proba[0][re2] >= 0)) {
      GraphUtils.linkElements(re1, re2, ranks, fathers);
      if (proba[0][re2] >= 0 && proba[0][re1] < 0) {
        for (k = 0; k < numLabels - 1; k++) {
          proba[k][re1] = proba[k][re2];
        }
      } else if (proba[0][re1] >= 0 && proba[0][re2] < 0) {
        for ( k = 0; k < numLabels - 1; k++) {
          proba[k][re2] = proba[k][re1];
        }
      }
    }
  }

  public Sequence getSegmentation() {
    return seqResult;
  }

  public Collection<ROI3DArea> getSegmentationROIs() {
    List<ROI3DArea> rois = new ArrayList<ROI3DArea>(colorLabels.length);
    for (int i = 0; i < colorLabels.length; i++) {
      rois.add(new ROI3DArea());
    }

    int[][] seqResData = seqResult.getDataXYZAsInt(0, 0);
    for (int z = 0; z < sz; z++) {
      for (int y = 0; y < sy; y++) {
        for (int x = 0; x < sx; x++) {
          rois.get((seqResData[z][x + y*sx] + 1) % colorLabels.length).addPoint(x, y, z);
        }
      }
    }

    for (int i = 0; i < colorLabels.length; i++) {
      rois.get(i).setColor(colorLabels[i]);
      rois.get(i).setName(nameLabels[i]);
    }

    return rois;
  }


}
