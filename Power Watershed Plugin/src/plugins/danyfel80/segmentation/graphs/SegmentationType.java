/**
 * 
 */
package plugins.danyfel80.segmentation.graphs;

/**
 * @author Daniel Felipe Gonzalez Obando Specifies the different algorithms used
 *         to make the segmentations.
 */
public enum SegmentationType {
  CollapseToSeeds("Collapse to seeds (p = q = 1)"), l2NormVoronoi(
      "l2 Norm Voronoi (p = 0, q = 2)"), l1NormVoronoi(
          "l1 Norm Voronoi (p = [0, 1], q = inf)"), GraphCuts(
              "Graph cuts (p = 1, q = 1)"), RandomWalker(
                  "Random walker (p = 1, q = 2)"), PowerWatershedQ1(
                      "Power watershed (p = inf, q = 1)"), PowerWatershedQ2(
                          "Power watershed (p = inf, q = 2)"), ShortestPathForest(
                              "Shortest path forest (p = inf, q = inf)");

  private final String name;

  SegmentationType(String name) {
    this.name = name;
  }

  @Override
  public String toString() {
    return name;
  }

  public static SegmentationType getAlgorithm(int p, int q) {
    if (p == 0) {
      if (q == 1) {
        return CollapseToSeeds;
      } else if (q == 2) {
        return l2NormVoronoi;
      } else {
        return l1NormVoronoi;
      }
    } else if (p > 0) {
      if (q == 1) {
        return GraphCuts;
      } else if (q == 2) {
        return RandomWalker;
      } else {
        return l1NormVoronoi;
      }
    } else {
      if (q == 1) {
        return PowerWatershedQ1;
      } else if (q == 2) {
        return PowerWatershedQ2;
      } else {
        return ShortestPathForest;
      }
    }
  }
}
