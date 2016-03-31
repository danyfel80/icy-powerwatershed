/**
 * 
 */
package algorithms.danyfel80.segmentation.graphcut;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class Edge {
  
  private Node   head;
  private Edge   nextEdge;
  private Edge   reverseEdge;
  private double residualCapacity;

  
  public Node getHead() {
    return head;
  }
  public void setHead(Node head) {
    this.head = head;
  }
  public Edge getNextEdge() {
    return nextEdge;
  }
  public void setNextEdge(Edge nextEdge) {
    this.nextEdge = nextEdge;
  }
  public Edge getReverseEdge() {
    return reverseEdge;
  }
  public void setReverseEdge(Edge reverseEdge) {
    this.reverseEdge = reverseEdge;
  }
  public double getResidualCapacity() {
    return residualCapacity;
  }
  public void setResidualCapacity(double residualCapacity) {
    this.residualCapacity = residualCapacity;
  }
  public void addResidualCapacity(double capacity) {
    this.residualCapacity += capacity;
  } 
}
