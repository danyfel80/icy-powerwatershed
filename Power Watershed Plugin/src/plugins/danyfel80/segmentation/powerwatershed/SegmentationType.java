/**
 * 
 */
package plugins.danyfel80.segmentation.powerwatershed;

/**
 * @author Daniel Felipe Gonzalez Obando
 * Specifies the different algorithms used to make the segmentations.
 */
public enum SegmentationType {
  CollapseToSeeds ("Collapse to seeds"),
  l2NormVoronoi ("l2 Norm Voronoi"),
  l1NormVoronoi ("l1 Norm Voronoi"),
  GraphCuts ("Graph cuts"),
  RandomWalker ("Random walker"),
  PowerWatershedQ1 ("Power watershed q = 1"),
  PowerWatershedQ2 ("Power watershed q = 2"),
  ShortestPathForest ("Shortest path forest");

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
