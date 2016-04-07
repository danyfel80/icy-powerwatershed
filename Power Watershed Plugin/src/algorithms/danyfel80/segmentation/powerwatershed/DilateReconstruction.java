package algorithms.danyfel80.segmentation.powerwatershed;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class DilateReconstruction {

  private static int MAX_W = 255;
  /**
   * Reconstruction by dilation of normalWeights under seedsFunction.  Union-find method described by Luc Vicent.
   * @param seedsFunction Image seeds
   * @param originalWeights Image weights
   * @param resultWeights Result of the reconstruction by dilation of normalWeights under seedsFunction.
   * @param edges List of couples of vertices forming edges
   * @param sx Size in X of sequence
   * @param sy Size in Y of sequence
   * @param sz Size in Z of sequence
   * @param M Amount of edges
   */
  public static void reconstruct(int[] seedsFunction, int[] originalWeights, 
      int[] resultWeights, int[][] edges, int sx, int sy,int sz, int M, boolean useQuickSort) {

    int k, p, i, n;
    boolean[] markedEdges = new boolean[M];
    int[] parentEdges = new int[M];
    int[] sortedEdges = new int[M];

    for (k = 0; k < M; k++) {
      parentEdges[k] = k;
      resultWeights[k] = seedsFunction[k];
      seedsFunction[k] = originalWeights[k];
      sortedEdges[k] = k;
    }

    if (useQuickSort) {
      GraphUtils.BinSort(seedsFunction, sortedEdges, M, MAX_W+1);
    } else {
      GraphUtils.quickStochasticSort(seedsFunction, sortedEdges, 0, M-1);
    }

    // First pass
    if (sz == 1) { // 2D image
      for (k = M - 1; k >= 0; k--) {
        p = sortedEdges[k];
        for (i = 1; i <= 6; i += 1) { // parcourt les 6 voisins  
          n = GraphUtils.getNeighborEdge(p, i, sx, sy, sz);
          if (n != -1) {
            if (markedEdges[n]) { 
              elementLinkGeodesicDilate(n, p, parentEdges, originalWeights, resultWeights);
            }
          }
          markedEdges[p]=true;
        }
      }
    } else { // 3D image
      for (k = M - 1; k >= 0; k--) {
        p = sortedEdges[k];
        for (i = 1; i <= 12; i += 1) { // parcourt les 12 voisins  
          n = GraphUtils.getNeighborEdge3D(edges[0][p], edges[1][p], p, i, sx, sy, sz);
          if (n != -1)
            if(markedEdges[n])
              elementLinkGeodesicDilate(n,p, parentEdges, originalWeights, resultWeights);
          markedEdges[p]=true;
        }
      }
    }

    // Second pass
    for(k = 0; k < M; k++) {
      p = sortedEdges[k];
      if (parentEdges[p] == p) { // p is root
        if (resultWeights[p] == MAX_W) {
          resultWeights[p] = originalWeights[p];
        }
      }
      else {
        resultWeights[p] = resultWeights[parentEdges[p]];
      }
    }
  }

  /**
   * 
   * @param childedge
   * @param parentEdge
   * @param fathers
   * @param normalWeights
   * @param weights
   */
  private static void elementLinkGeodesicDilate(int childedge, int parentEdge, int[] fathers,
      int[] normalWeights, int[] weights) {
    int reconstructedNode = GraphUtils.findElement(childedge, fathers);

    if (reconstructedNode != parentEdge) {
      if(normalWeights[reconstructedNode] == normalWeights[parentEdge] || normalWeights[parentEdge] >= weights[reconstructedNode]) {
        fathers[reconstructedNode] = parentEdge;
        weights[parentEdge] = Math.max(weights[reconstructedNode], weights[parentEdge]);
      } else {
        weights[parentEdge] = MAX_W;
      }
    }
  }
}
