/**
 * 
 */
package algorithms.danyfel80.segmentation.graphcut;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class Node {
  private int    id;
  private Edge firstOutgoingEdge;
  private Edge parentEdge;
  private Node nextActiveNode;
  
  private long      timeStamp;
  private int       distanceToTerminalNode;
  private boolean   inSinkTree;
  private double    terminalResidualCapacity;
  
  
  public int getId() {
    return id;
  }
  public void setId(int id) {
    this.id = id;
  }
  
  public Edge getFirstOutgoingEdge() {
    return firstOutgoingEdge;
  }
  public void setFirstOutgoingEdge(Edge firstOutgoingEdge) {
    this.firstOutgoingEdge = firstOutgoingEdge;
  }
  public Edge getParentEdge() {
    return parentEdge;
  }
  public void setParentEdge(Edge parentEdge) {
    this.parentEdge = parentEdge;
  }
  public Node getNextActiveNode() {
    return nextActiveNode;
  }
  public void setNextActiveNode(Node nextActiveNode) {
    this.nextActiveNode = nextActiveNode;
  }
  public long getTimeStamp() {
    return timeStamp;
  }
  public void setTimeStamp(long timeStamp) {
    this.timeStamp = timeStamp;
  }
  public int getDistanceToTerminalNode() {
    return distanceToTerminalNode;
  }
  public void setDistanceToTerminalNode(int distanceToTerminalNode) {
    this.distanceToTerminalNode = distanceToTerminalNode;
  }
  public boolean isInSinkTree() {
    return inSinkTree;
  }
  public void setInSinkTree(boolean belongsToSink) {
    this.inSinkTree = belongsToSink;
  }
  public double getTerminalResidualCapacity() {
    return terminalResidualCapacity;
  }
  public void setTerminalResidualCapacity(double terminalResidualCapacity) {
    this.terminalResidualCapacity = terminalResidualCapacity;
  }
  public void addTerminalResidualCapacity(double capacity) {
    this.terminalResidualCapacity += capacity;
  }
  public boolean isActive() {
    return this.nextActiveNode != null;
  }
  
}
