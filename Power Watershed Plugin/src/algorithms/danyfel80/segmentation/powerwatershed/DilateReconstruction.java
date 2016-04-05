package algorithms.danyfel80.segmentation.powerwatershed;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class DilateReconstruction {

  private static int MAX_W;
  /**
   * Reconstruction by dilation of normalWeights under seedsFunction.  Union-find method described by Luc Vicent.
   */
  public static int[] reconstruct(int[] seedsFunction, int[] normalWeights, 
      int[] weights, int[][] edges, int sx, int sy,int sz, int M) {

    int k, p, i, n;
    boolean[] marked = new boolean[M];
    int[] fathers = new int[M];
    int[] edgesSorted = new int[M];

    for (k = 0; k < M; k++) {
      fathers[k] = k;
      weights[k] = seedsFunction[k];
      seedsFunction[k] = normalWeights[k];
      edgesSorted[k] = k;
    }

    MAX_W = 255;
    GraphUtils.sort(seedsFunction, edgesSorted, M, MAX_W+1);

    // First pass
    if (sz == 1) { // 2D image
      for (k = M - 1; k >= 0; k--) {
        p = edgesSorted[k];
        for (i = 1; i <= 6; i += 1) { // parcourt les 6 voisins  
          n = GraphUtils.getNeighborEdge(p, i, sx, sy, sz);
          if (n != -1) {
            if (marked[n]) { 
              elementLinkGeodesicDilate(n,p, fathers, normalWeights, weights);
            }
          }
          marked[p]=true;
        }
      }
    } else { // 3D image
      for (k = M - 1; k >= 0; k--) {
        p = edgesSorted[k];
        for (i = 1; i <= 12; i += 1) { // parcourt les 12 voisins  
          n = GraphUtils.getNeighborEdge3D(edges[0][p], edges[1][p], p, i, sx, sy, sz);
          if (n != -1)
            if(marked[n])
              elementLinkGeodesicDilate(n,p, fathers, normalWeights, weights);
          marked[p]=true;
        }
      }
    }

    // Second pass
    for(k = 0; k < M; k++) {
      p = edgesSorted[k];
      if (fathers[p] == p) { // p is root
        if (weights[p] == MAX_W) {
          weights[p] = normalWeights[p];
        }
      }
      else {
        normalWeights[p] = normalWeights[fathers[p]];
      }
    }
    return weights;
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
