package plugins.danyfel80.segmentation.powerwatershed.classes.graphcut;

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
  private List<Edge> edges;

  private double flow = 0;

  private LinkedList<Node>[] activeNodes;
  private LinkedList<Node> orphanNodes;
  
  private long TIME;



  public Graph(int numNodes, int numEdges) {
    nodes = new ArrayList<Node>(numNodes);
    edges = new ArrayList<Edge>(2*numEdges);
  }

  public int addNodes(int amount) {
    int nodePos = nodes.size();
    for (int i = 0; i < amount; i++)
      nodes.add(new Node());
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
    
    System.out.println("initialized");
    
    Node p, q, currentNode;
    p = q = currentNode = null;
    Edge e = null;


    while (true) {

      // Test consistency
      p = currentNode;
      if (p != null) {
        p.setNextActiveNode(null);
        if(p.getParentEdge() == null) {
          p = null;
        }
      }

      if (p == null) {
        p = getNextActiveNode();
        if (p == null) {
          break;
        }
      }

      // growTrees
      if(!p.isInSinkTree()) {
        System.out.println("grow source");
        // grow source tree
        for (e = p.getFirstOutgoingEdge(); e != null; e=e.getNextEdge()) {
          if (e.getResidualCapacity() != 0.0) {
            q = e.getHead();
            if (q.getParentEdge() == null) {
              q.setInSinkTree(false);
              q.setParentEdge(e.getReverseEdge());
              q.setTimeStamp(p.getTimeStamp());
              q.setDistanceToTerminalNode(q.getDistanceToTerminalNode()+1);
              setActiveNode(q);
            }
            else if (q.isInSinkTree()) {
              break;
            }
            else if (q.getTimeStamp() <= p.getTimeStamp() &&
                q.getDistanceToTerminalNode() > p.getDistanceToTerminalNode()) {
              q.setParentEdge(e.getReverseEdge());
              q.setTimeStamp(p.getTimeStamp());
              q.setDistanceToTerminalNode(p.getDistanceToTerminalNode()+1);
            }
          }
        }
      }
      else {
        System.out.println("grow sink");
        // grow sink tree
        for (e = p.getFirstOutgoingEdge(); e != null; e=e.getNextEdge()) {
          if (e.getReverseEdge().getResidualCapacity() != 0.0) {
            q = e.getHead();
            if (q.getParentEdge() == null) {
              q.setInSinkTree(true);
              q.setParentEdge(e.getReverseEdge());
              q.setTimeStamp(p.getTimeStamp());
              q.setDistanceToTerminalNode(q.getDistanceToTerminalNode()+1);
              setActiveNode(q);
            }
            else if (!q.isInSinkTree()) {
              e = e.getReverseEdge();
              break;
            }
            else if (q.getTimeStamp() <= p.getTimeStamp() && 
                q.getDistanceToTerminalNode() > p.getDistanceToTerminalNode()) {
              q.setParentEdge(e.getReverseEdge());
              q.setTimeStamp(p.getTimeStamp());
              q.setDistanceToTerminalNode(p.getDistanceToTerminalNode()+1);
            }
          }
        }
      }
    }
    
    TIME++;
    
    System.out.println("growing phase stopped");
    
    if (e != null) {
      System.out.println("augmenting path found");
      // augmenting path found
      p.setNextActiveNode(p);
      currentNode = p;

      // Augment flow
      augmentFlow(e);

      System.out.println("flow saturated");
      
      // Adopt orphans
      while (!orphanNodes.isEmpty()) {
        Node o = orphanNodes.element();
        o.setNextActiveNode(null);

        while (!orphanNodes.isEmpty()) {
          o = orphanNodes.element();
          orphanNodes.poll();
          p = o;
          if (p.isInSinkTree())
            adoptSinkOrphan(p);
          else
            adoptSourceOrphan(p);
        }
        orphanNodes.poll();
      }
      
      System.out.println("orphans adopted");
    }
    else {
      currentNode = null;
    }

    return 0;
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
      if (!activeNodes[1].isEmpty()) {
        activeNodes[1].getLast().setNextActiveNode(node);
      }
      node.setNextActiveNode(node);
      activeNodes[1].add(node);
    }
  }

  private Node getNextActiveNode() {
    Node n;

    while (true) {
      n = activeNodes[0].peek();
      if (n == null) {
        LinkedList<Node> temp = activeNodes[0];
        activeNodes[0] = activeNodes[1];
        activeNodes[1] = temp;
        n = activeNodes[0].peek();
        if (n == null)
          return null;
      }

      if (n.getNextActiveNode() == n) {
        activeNodes[0].remove();
      }
      n.setNextActiveNode(null);

      if (n.getParentEdge() != null) {
        return n;
      }
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

    // 2.a. Source tree
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
    Edge e0, e0Min = null, e;
    int d, dMin = Integer.MAX_VALUE;

    // Try to find a new parent
    for (e0 = p.getFirstOutgoingEdge(); e0 != null; e0 = e0.getNextEdge()) {
      if (e0.getReverseEdge().getResidualCapacity() > 0.0) {
        q = e0.getHead();
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
                e0Min = e0;
                dMin = d;
              }
              // set marks along the path
              for (q = e0.getHead(); q.getTimeStamp() != TIME; q = q.getParentEdge().getHead()) {
                q.setTimeStamp(TIME);
                q.setDistanceToTerminalNode(d--);
              }
              
            }
          }
        }
      }
    }
    
    p.setParentEdge(e0Min);
    if (p.getParentEdge() != null) {
      p.setTimeStamp(TIME);
      p.setDistanceToTerminalNode(dMin + 1);
    }
    else {
      // No parent was found
      // Process neighbors
      for(e0 = p.getFirstOutgoingEdge(); e0 != null; e0 = e0.getNextEdge()) {
        q = e0.getHead();
        
        if (!q.isInSinkTree()) {
          e = q.getParentEdge();
          if (e != null) {
            if (e0.getReverseEdge().getResidualCapacity() != 0.0)
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
    Edge e0, e0Min = null, e;
    int d, dMin = Integer.MAX_VALUE;

    // Try to find a new parent
    for (e0 = p.getFirstOutgoingEdge(); e0 != null; e0 = e0.getNextEdge()) {
      if (e0.getResidualCapacity() > 0.0) {
        q = e0.getHead();
        
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
                e0Min = e0;
                dMin = d;
              }
              // set marks along the path
              for (q = e0.getHead(); q.getTimeStamp() != TIME; q = q.getParentEdge().getHead()) {
                q.setTimeStamp(TIME);
                q.setDistanceToTerminalNode(d--);
              }
              
            }
          }
        }
      }
    }
    
    p.setParentEdge(e0Min);
    if (p.getParentEdge() != null) {
      p.setTimeStamp(TIME);
      p.setDistanceToTerminalNode(dMin + 1);
    }
    else {
      for(e0 = p.getFirstOutgoingEdge(); e0 != null; e0 = e0.getNextEdge()) {
        q = e0.getHead();
        
        if (!q.isInSinkTree()) {
          e = q.getParentEdge();
          if (e != null) {
            if (e0.getReverseEdge().getResidualCapacity() != 0.0)
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



  @SuppressWarnings("unchecked")
  private void initMaxFlow() {
    activeNodes = new LinkedList[2];
    activeNodes[0] = new LinkedList<Node>();
    activeNodes[1] = new LinkedList<Node>();
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
