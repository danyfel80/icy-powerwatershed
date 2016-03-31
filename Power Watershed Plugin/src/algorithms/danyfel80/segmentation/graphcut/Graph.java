package algorithms.danyfel80.segmentation.graphcut;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Daniel Felipe Gonzalez Obando
 *
 */
public class Graph {

  public enum TerminalType {
    SOURCE,
    SINK
  }

  private static Edge TERMINAL = new Edge();
  private static Edge ORPHAN = new Edge();

  private List<Node> nodes;
  public List<Edge> edges;

  private double flow = 0;

  private LinkedList<Node> activeNodes;
  private LinkedList<Node> orphanNodes;
  
  private long TIME;



  public Graph(int numNodes, int numEdges) {
    nodes = new ArrayList<Node>(numNodes);
    edges = new ArrayList<Edge>(2*numEdges);
  }

  public int addNodes(int amount) {
    int nodePos = nodes.size();
    for (int i = 0; i < amount; i++){
      Node n = new Node();
      n.setId(nodePos+i);
      nodes.add(n);
      
    }
    return nodePos;
  }

  public void addEdge(int nodePId, int nodeQId, double capacityPQ, double capacityQP) {
    assert(nodePId >= 0 && nodePId < this.nodes.size());
    assert(nodeQId >= 0 && nodeQId < this.nodes.size());
    assert(nodePId != nodeQId);
    assert(capacityPQ >= 0.0);
    assert(capacityQP >= 0.0);

    Edge e = new Edge();
    Edge eR = new Edge();
    edges.add(e);
    edges.add(eR);

    Node p = nodes.get(nodePId);
    Node q = nodes.get(nodeQId);

    e.setReverseEdge(eR);
    e.setNextEdge(p.getFirstOutgoingEdge());
    e.setHead(q);
    e.setResidualCapacity(capacityPQ);
    p.setFirstOutgoingEdge(e);

    eR.setReverseEdge(e);
    eR.setNextEdge(q.getFirstOutgoingEdge());
    eR.setHead(p);
    eR.setResidualCapacity(capacityQP);
    q.setFirstOutgoingEdge(eR);
  }

  public void addTerminalWeights(int nodeId, double capacitySource, double capacitySink) {
    assert(nodeId >= 0 && nodeId < this.nodes.size());

    Node node = nodes.get(nodeId);

    double delta = node.getTerminalResidualCapacity();
    if (delta > 0)
      capacitySource += delta;
    else
      capacitySink -= delta;

    flow += ((capacitySource < capacitySink)? capacitySource: capacitySink);
    node.setTerminalResidualCapacity(capacitySource - capacitySink);
  }

  public double computeMaxFlow() {
    initMaxFlow();
    
    //System.out.println("initialized");
    
    Node activeNode, headNode, currentNode;
    activeNode = headNode = currentNode = null;
    Edge e = null;


    while (true) {

      activeNode = currentNode;
      if (activeNode != null) {
        activeNode.setNextActiveNode(null);
        if(activeNode.getParentEdge() == null) {
          activeNode = null;
        }
      }

      if (activeNode == null) {
        activeNode = getNextActiveNode();
        if (activeNode == null) {
          //System.out.println("no next active node");
          break;
        }
      }
      
      //System.out.println("on node " + activeNode.getId());
      // growTrees
      if(!activeNode.isInSinkTree()) {
        //System.out.println("grow source");
        // grow source tree
        for (e = activeNode.getFirstOutgoingEdge(); e != null; e=e.getNextEdge()) {
          if (e.getResidualCapacity() != 0.0) {
            headNode = e.getHead();
            
            if (headNode.getParentEdge() == null) {
              headNode.setInSinkTree(false);
              headNode.setParentEdge(e.getReverseEdge());
              headNode.setTimeStamp(activeNode.getTimeStamp());
              headNode.setDistanceToTerminalNode(activeNode.getDistanceToTerminalNode()+1);
              setActiveNode(headNode);
            }
            else if (headNode.isInSinkTree()) {
              break;
            }
            else if (headNode.getTimeStamp() <= activeNode.getTimeStamp() &&
                headNode.getDistanceToTerminalNode() > activeNode.getDistanceToTerminalNode()) {
              headNode.setParentEdge(e.getReverseEdge());
              headNode.setTimeStamp(activeNode.getTimeStamp());
              headNode.setDistanceToTerminalNode(activeNode.getDistanceToTerminalNode()+1);
            }
          }
        }
      }
      else {
        //System.out.println("grow sink");
        // grow sink tree
        for (e = activeNode.getFirstOutgoingEdge(); e != null; e=e.getNextEdge()) {
          if (e.getReverseEdge().getResidualCapacity() != 0.0) {
            headNode = e.getHead();
            if (headNode.getParentEdge() == null) {
              headNode.setInSinkTree(true);
              headNode.setParentEdge(e.getReverseEdge());
              headNode.setTimeStamp(activeNode.getTimeStamp());
              headNode.setDistanceToTerminalNode(activeNode.getDistanceToTerminalNode()+1);
              setActiveNode(headNode);
            }
            else if (!headNode.isInSinkTree()) {
              e = e.getReverseEdge();
              break;
            }
            else if (headNode.getTimeStamp() <= activeNode.getTimeStamp() && 
                headNode.getDistanceToTerminalNode() > activeNode.getDistanceToTerminalNode()) {
              headNode.setParentEdge(e.getReverseEdge());
              headNode.setTimeStamp(activeNode.getTimeStamp());
              headNode.setDistanceToTerminalNode(activeNode.getDistanceToTerminalNode()+1);
            }
          }
        }
      }
    }
    
    TIME++;
    
    //System.out.println("growing phase stopped");
    
    if (e != null) {
      //System.out.println("augmenting path found");
      // augmenting path found
      if (activeNode != null)
        activeNode.setNextActiveNode(activeNode);
      currentNode = activeNode;

      // Augment flow
      augmentFlow(e);

      //System.out.println("flow saturated");
      
      // Adopt orphans
      while (!orphanNodes.isEmpty()) {
        Node o = orphanNodes.poll();
        o.setNextActiveNode(null);
        if (o.isInSinkTree())
          adoptSinkOrphan(o);
        else
          adoptSourceOrphan(o);
      }
      
      //System.out.println("orphans adopted");
    }
    else {
      currentNode = null;
    }

    return flow;
  }

  public TerminalType getSegment(int nodeId) {
    assert(nodeId >= 0 && nodeId < this.nodes.size());
    if (nodes.get(nodeId).getParentEdge() != null) {
      return nodes.get(nodeId).isInSinkTree()? TerminalType.SINK: TerminalType.SOURCE;
    }
    else {
      return null;
    }
  }


  private void setActiveNode(Node node) {
    if (!node.isActive()) {
      /*if (!activeNodes[1].isEmpty()) {
        activeNodes[1].getLast().setNextActiveNode(node);
      }
      node.setNextActiveNode(node);
      activeNodes[1].add(node);*/
      if (!activeNodes.isEmpty()) {
        activeNodes.getLast().setNextActiveNode(node);
      }
      node.setNextActiveNode(node);
      activeNodes.add(node);
    }
  }

  private Node getNextActiveNode() {
    Node n;

    while (true) {
      n = activeNodes.poll();
      if (n == null) {
        return null;
      }
  
      n.setNextActiveNode(null);
  
      if (n.getParentEdge() != null) {
        return n;
      }
//      n = activeNodes[0].peek();
//      if (n == null) {
//        LinkedList<Node> temp = activeNodes[0];
//        activeNodes[0] = activeNodes[1];
//        activeNodes[1] = temp;
//        n = activeNodes[0].peek();
//        if (n == null)
//          return null;
//      }
//
//      if (n.getNextActiveNode() == n) {
//        activeNodes[0].clear();
//      }
//      else {
//        activeNodes[0].remove();
//      }
//      n.setNextActiveNode(null);
//
//      if (n.getParentEdge() != null) {
//        return n;
//      }
    }
  }

  private void augmentFlow(Edge middleE) {
    Node p;
    Edge e = null;
    double bottleneck;

    // 1. Find bottleneck capacity
    // 1.a. Source tree
    bottleneck = middleE.getResidualCapacity();
    for (p = middleE.getReverseEdge().getHead(); ; p = e.getHead()) {
      e = p.getParentEdge();
      if (e == TERMINAL) break;
      if (bottleneck > e.getReverseEdge().getResidualCapacity()) {
        bottleneck = e.getReverseEdge().getResidualCapacity();
      }
    }

    if (bottleneck > p.getTerminalResidualCapacity())
      bottleneck = p.getTerminalResidualCapacity();

    // 1.b. Sink tree
    for (p = middleE.getHead(); ; p = e.getHead()) {
      e = p.getParentEdge();
      if (e == TERMINAL) break;
      if (bottleneck > e.getResidualCapacity()) {
        bottleneck = e.getResidualCapacity();
      }
    }

    if (bottleneck > -p.getTerminalResidualCapacity())
      bottleneck = -p.getTerminalResidualCapacity();

    // 2. Augment flow
    middleE.getReverseEdge().addResidualCapacity(bottleneck);
    middleE.addResidualCapacity(-bottleneck);
    // 2.a. Source tree
    for (p = middleE.getReverseEdge().getHead(); ; p = e.getHead()) {
      e = p.getParentEdge();
      if (e == TERMINAL) 
        break;

      e.addResidualCapacity(bottleneck);
      e.getReverseEdge().addResidualCapacity(-bottleneck);
      if(e.getReverseEdge().getResidualCapacity() == 0.0) {
        addFrontOrphan(p);
      }
    }
    p.addTerminalResidualCapacity(-bottleneck);
    if (p.getTerminalResidualCapacity() == 0.0) {
      addFrontOrphan(p);
    }

    // 2.b. Sink tree
    for (p = middleE.getHead(); ; p = e.getHead()) {
      e = p.getParentEdge();
      if (e == TERMINAL) 
        break;

      e.addResidualCapacity(-bottleneck);
      e.getReverseEdge().addResidualCapacity(bottleneck);
      if(e.getResidualCapacity() == 0.0) {
        addFrontOrphan(p);
      }
    }
    p.addTerminalResidualCapacity(bottleneck);
    if (p.getTerminalResidualCapacity() == 0.0) {
      addFrontOrphan(p);
    }

    flow += bottleneck;
  }

  private void adoptSourceOrphan(Node p) {
    Node q;
    Edge orphanEdge, orphanEdgeMin = null, e;
    int d, dMin = Integer.MAX_VALUE;

    // Try to find a new parent
    for (orphanEdge = p.getFirstOutgoingEdge(); orphanEdge != null; orphanEdge = orphanEdge.getNextEdge()) {
      if (orphanEdge.getReverseEdge().getResidualCapacity() != 0.0) {
        q = orphanEdge.getHead();
        if (!q.isInSinkTree()) {
          e = q.getParentEdge();
          if(e != null) {
            // check origin
            d = 0;
            while (true) {
              if (q.getTimeStamp() == TIME) {
                d += q.getDistanceToTerminalNode();
                break;
              }
              e = q.getParentEdge();
              d++;
              if (e == TERMINAL) {
                q.setTimeStamp(TIME);
                q.setDistanceToTerminalNode(1);
                break;
              }
              if (e == ORPHAN) {
                d = Integer.MAX_VALUE;
                break;
              }
              q = e.getHead();
            }
            
            if (d < Integer.MAX_VALUE) { // q originates from the source
              if (d < dMin) {
                orphanEdgeMin = orphanEdge;
                dMin = d;
              }
              // set marks along the path
              for (q = orphanEdge.getHead(); q.getTimeStamp() != TIME; q = q.getParentEdge().getHead()) {
                q.setTimeStamp(TIME);
                q.setDistanceToTerminalNode(d);
                d--;
              }
              
            }
          }
        }
      }
    }
    
    p.setParentEdge(orphanEdgeMin);
    if (p.getParentEdge() != null) {
      p.setTimeStamp(TIME);
      p.setDistanceToTerminalNode(dMin + 1);
    }
    else {
      // No parent was found
      // Process neighbors
      for(orphanEdge = p.getFirstOutgoingEdge(); orphanEdge != null; orphanEdge = orphanEdge.getNextEdge()) {
        q = orphanEdge.getHead();
        
        if (!q.isInSinkTree()) {
          e = q.getParentEdge();
          if (e != null) {
            if (orphanEdge.getReverseEdge().getResidualCapacity() != 0.0)
              setActiveNode(q);
            if (e != TERMINAL && e != ORPHAN && e.getHead() == p) {
              addBackOrphan(q); // add q to end of orphans list
            }
          }
        }
      }
    }
  }

  private void adoptSinkOrphan(Node p) {
    Node q;
    Edge orphanEdge, orphanEdgeMin = null, e;
    int d, dMin = Integer.MAX_VALUE;

    // Try to find a new parent
    for (orphanEdge = p.getFirstOutgoingEdge(); orphanEdge != null; orphanEdge = orphanEdge.getNextEdge()) {
      if (orphanEdge.getResidualCapacity() != 0.0) {
        q = orphanEdge.getHead();
        
        if (!q.isInSinkTree()) {
          e = q.getParentEdge();
          if(e != null) {
            // check origin
            d = 0;
            while (true) {
              if (q.getTimeStamp() == TIME) {
                d += q.getDistanceToTerminalNode();
                break;
              }
              e = q.getParentEdge();
              d++;
              if (e == TERMINAL) {
                q.setTimeStamp(TIME);
                q.setDistanceToTerminalNode(1);
                break;
              }
              if (e == ORPHAN) {
                d = Integer.MAX_VALUE;
                break;
              }
              q = e.getHead();
            }
            
            if (d < Integer.MAX_VALUE) { // q originates from the source
              if (d < dMin) {
                orphanEdgeMin = orphanEdge;
                dMin = d;
              }
              // set marks along the path
              for (q = orphanEdge.getHead(); q.getTimeStamp() != TIME; q = q.getParentEdge().getHead()) {
                q.setTimeStamp(TIME);
                q.setDistanceToTerminalNode(d);
                d--;
              }
              
            }
          }
        }
      }
    }
    
    p.setParentEdge(orphanEdgeMin);
    if (p.getParentEdge() != null) {
      p.setTimeStamp(TIME);
      p.setDistanceToTerminalNode(dMin + 1);
    }
    else {
      for(orphanEdge = p.getFirstOutgoingEdge(); orphanEdge != null; orphanEdge = orphanEdge.getNextEdge()) {
        q = orphanEdge.getHead();
        
        if (!q.isInSinkTree()) {
          e = q.getParentEdge();
          if (e != null) {
            if (orphanEdge.getReverseEdge().getResidualCapacity() != 0.0)
              setActiveNode(q);
            if (e != TERMINAL && e != ORPHAN && e.getHead() == p) {
              addBackOrphan(q);
            }
          }
        }
      }
    }
  }

  private void addFrontOrphan(Node n) {
    n.setParentEdge(ORPHAN);
    n.setNextActiveNode(orphanNodes.peek());
    orphanNodes.addFirst(n);
  }

  private void addBackOrphan(Node n) {
    n.setParentEdge(ORPHAN);
    if (!orphanNodes.isEmpty()) orphanNodes.getLast().setNextActiveNode(n);
    orphanNodes.addLast(n);
  }



  private void initMaxFlow() {
    activeNodes = new LinkedList<Node>();
    //activeNodes = new LinkedList[2];
    //activeNodes[0] = new LinkedList<Node>();
    //activeNodes[1] = new LinkedList<Node>();
    orphanNodes = new LinkedList<Node>();

    TIME = 0;
    
    for (Node n : nodes) {
      n.setNextActiveNode(null);
      n.setTimeStamp(TIME);
      if(n.getTerminalResidualCapacity() > 0) {
        n.setInSinkTree(false);
        n.setParentEdge(TERMINAL);
        setActiveNode(n);
        n.setDistanceToTerminalNode(1);
      }
      else if (n.getTerminalResidualCapacity() < 0) {
        n.setInSinkTree(true);
        n.setParentEdge(TERMINAL);
        setActiveNode(n);
        n.setDistanceToTerminalNode(1);
      }
      else {
        n.setParentEdge(null);
      }
    }
  }

  public double getFlow() {
    return flow;
  }

  public void reset() {
    nodes.clear();
    edges.clear();
    activeNodes = null;
    orphanNodes = null;
    flow = 0;
    TIME = 0;
    System.gc();
  }
}
