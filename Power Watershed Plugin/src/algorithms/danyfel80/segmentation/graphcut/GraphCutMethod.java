package algorithms.danyfel80.segmentation.graphcut;

import java.util.LinkedList;
import java.util.List;

/**
 * Class implementing the grach cut algorithm.
 */
public class GraphCutMethod {

  // graph structure
  private GraphGC graph;

  // counter for initialisation of edges
  private int edgeNum;

  // the total flow in the whole graph
  private float totalFlow;

  // counter for the numbers of iterations to maxflow
  private int maxflowIteration;

  // Lists of active nodes: activeQueueFirst points to first
  // elements of the lists, activeQueueLast to the last ones.
  // In between, nodes are connected via reference to next node
  // in each node.
  private int[] activeQueueFirst;
  private int[] activeQueueLast;

  // list of orphans
  private LinkedList<Integer> orphans;

  // counter for iterations of main loop
  private int time;

  public int maxTerminalWeight = 50;
  
  /**
   * Initialises the graph cut implementation and allocates the memory needed
   * for the given number of nodes and edges.
   *
   * @param numNodes The number of nodes that should be created.
   * @param numEdges The number of edges that you can add. A directed edge and its 
   *                 counterpart (i.e., the directed edge in the other
   *                 direction) count as one edge.
   */
  public GraphCutMethod(int numNodes, int numEdges) {
    graph            = new GraphGC(numNodes, numEdges);
    edgeNum          = 0;
    totalFlow        = 0;
    maxflowIteration = 0;
    activeQueueFirst = new int[2];
    activeQueueLast  = new int[2];
    orphans          = new LinkedList<Integer>();
  }

  /**
   * Set the affinity for one node to belong to the foreground (i.e., source)
   * or background (i.e., sink).
   *
   * @param node   The number of the node.
   * @param source The affinity of this node to the foreground (i.e., source)
   * @param sink   The affinity of this node to the background (i.e., sink)
   */
  public void setTerminalWeights(int node, float source, float sink, boolean update) {
    
    source = (source > maxTerminalWeight)? maxTerminalWeight: source;
    sink = (sink > maxTerminalWeight)? maxTerminalWeight: sink;
    
    if(update) {
      float delta = graph.getResidualNodeCapacity(node);
    
      if (delta > 0)
        source += delta;
      else
        sink   -= delta;
    }
    
    
    totalFlow += (source < sink) ? source : sink;
    graph.setResidualNodeCapacity(node, source - sink);
  
  }

  /**
   * Set the edge weight of an undirected edge between two nodes.
   *
   * Please note that you cannot call any <tt>setEdgeWeight</tt> more often
   * than the number of edges you specified at the time of construction!
   *
   * @param node1   The first node.
   * @param node2   The second node.
   * @param weight  The weight (i.e., the cost) of the connecting edge.
   */
  public void setEdgeWeight(int node1, int node2, float weight) {

    setEdgeWeight(node1, node2, weight, weight);
  }

  /**
   * Set the edge weight of a pair of directed edges between two nodes.
   *
   * Please note that you cannot call any <tt>setEdgeWeight</tt> more often
   * than the number of edges you specified at the time of construction!
   *
   * @param node1      The first node.
   * @param node2      The second node.
   * @param weight1to2 The weight (i.e., the cost) of the directed edge from
   *                   node1 to node2.
   * @param weight2to1 The weight (i.e., the cost) of the directed edge from
   *                   node2 to node1.
   */
  public void setEdgeWeight(int node1, int node2, float weight1to2, float weight2to1) {

    // get edge indices
    int edge        = edgeNum; edgeNum++;
    int reverseEdge = edgeNum; edgeNum++;

    // link edges
    graph.setSister(edge, reverseEdge);
    graph.setSister(reverseEdge, edge);

    // add node1 to edge
    graph.setNextEdge(edge, graph.getFirstOutgoing(node1));
    graph.setFirstOutgoing(node1, edge);

    // add node2 to reverseEdge
    graph.setNextEdge(reverseEdge, graph.getFirstOutgoing(node2));
    graph.setFirstOutgoing(node2, reverseEdge);

    // set targets of edges
    graph.setHead(edge, node2);
    graph.setHead(reverseEdge, node1);

    // set residual capacities
    graph.setResidualEdgeCapacity(edge, weight1to2);
    graph.setResidualEdgeCapacity(reverseEdge, weight2to1);
  }

  /**
   * Performs the actual max-flow/min-cut computation.
   *
   * @param reuseTrees   reuse trees of a previos call
   * @param changedNodes list of nodes that potentially changed their
   *                     segmentation compared to a previous call, can be set
   *                     to <tt>null</tt>
   */
  public float computeMaximumFlow(boolean reuseTrees, List<Integer> changedNodes) {

    if (maxflowIteration == 0)
      reuseTrees = false;

    if (reuseTrees)
      maxflowReuseTreesInit();
    else
      maxflowInit();

    int currentNode = GraphGC.NONE;
    int edge        = GraphGC.NONE;

    // main loop
    while (true) {

      int activeNode = currentNode;

      if (activeNode != GraphGC.NONE) {
        // remove active flag
        graph.setNextNode(activeNode, GraphGC.NONE);
        if (graph.getParent(activeNode) == GraphGC.NONE)
          activeNode = GraphGC.NONE;
      }
      if (activeNode == GraphGC.NONE) {
        activeNode = getNextActiveNode();
        if (activeNode == GraphGC.NONE)
          // no more active nodes - we're done here
          break;
      }

      // groth
      if (!graph.isInSink(activeNode)) {
        // grow source tree
        for (edge = graph.getFirstOutgoing(activeNode); edge != GraphGC.NONE; edge = graph.getNextEdge(edge)) {
          if (graph.getResidualEdgeCapacity(edge) != 0) {

            int headNode = graph.getHead(edge);

            if (graph.getParent(headNode) == GraphGC.NONE) {
              // free node found, add to source tree
              graph.isInSink(headNode, false);
              graph.setParent(headNode, graph.getSister(edge));
              graph.setTimestamp(headNode, graph.getTimestamp(activeNode));
              graph.setDistance(headNode, graph.getDistance(activeNode) + 1);
              setNodeActive(headNode);
              addToChangedList(headNode);

            } else if (graph.isInSink(headNode)) {
              // node is not free and belongs to other tree - path
              // via edge found
              break;

            } else if (graph.getTimestamp(headNode) <= graph.getTimestamp(activeNode) &&
                graph.getDistance(headNode)  >  graph.getDistance(activeNode)) {
              // node is not free and belongs to our tree - try to
              // shorten its distance to the source
              graph.setParent(headNode, graph.getSister(edge));
              graph.setTimestamp(headNode, graph.getTimestamp(activeNode));
              graph.setDistance(headNode, graph.getDistance(activeNode) + 1);
            }
          }
        }
      } else {
        // activeNode is in sink, grow sink tree
        for (edge = graph.getFirstOutgoing(activeNode); edge != GraphGC.NONE; edge = graph.getNextEdge(edge)) {
          if (graph.getResidualEdgeCapacity(graph.getSister(edge)) != 0) {

            int headNode = graph.getHead(edge);

            if (graph.getParent(headNode) == GraphGC.NONE) {
              // free node found, add to sink tree
              graph.isInSink(headNode, true);
              graph.setParent(headNode, graph.getSister(edge));
              graph.setTimestamp(headNode, graph.getTimestamp(activeNode));
              graph.setDistance(headNode, graph.getDistance(activeNode) + 1);
              setNodeActive(headNode);
              addToChangedList(headNode);

            } else if (!graph.isInSink(headNode)) {
              // node is not free and belongs to other tree - path
              // via edge's sister found
              edge = graph.getSister(edge);
              break;

            } else if (graph.getTimestamp(headNode) <= graph.getTimestamp(activeNode) &&
                graph.getDistance(headNode)  >  graph.getDistance(activeNode)) {
              // node is not free and belongs to our tree - try to
              // shorten its distance to the sink
              graph.setParent(headNode, graph.getSister(edge));
              graph.setTimestamp(headNode, graph.getTimestamp(activeNode));
              graph.setDistance(headNode, graph.getDistance(activeNode) + 1);
            }
          }
        }
      }

      time++;

      if (edge != GraphGC.NONE) {
        // we found a path via edge

        // set active flag
        graph.setNextNode(activeNode, activeNode);
        currentNode = activeNode;

        // augmentation
        augment(edge);

        // adoption
        while (orphans.size() > 0) {
          int orphan = orphans.poll();
          if (graph.isInSink(orphan))
            processSinkOrphan(orphan);
          else
            processSourceOrphan(orphan);
        }
      } else {
        // no path found
        currentNode = GraphGC.NONE;
      }
    }

    maxflowIteration++;

    // create list of changed nodes
    if (changedNodes != null) {
      changedNodes.clear();
      for (int i = 0; i < graph.getNumNodes(); i++)
        if (graph.isInChangedList(i))
          changedNodes.add(i);
    }

    return totalFlow;
  }

  /**
   * Get the segmentation, i.e., the terminal node that is connected to the
   * specified node. If there are several min-cut solutions, free nodes are
   * assigned to the background.
   *
   * @param nodeId the node to check
   * @return Either <tt>Terminal.FOREGROUND</tt> or
   *         <tt>Terminal.BACKGROUND</tt>
   */
  public Terminal getTerminal(int node) {

    if (graph.getParent(node) != GraphGC.NONE)
      return graph.isInSink(node) ? Terminal.BACKGROUND : Terminal.FOREGROUND;
    else
      return Terminal.BACKGROUND;
  }

  /**
   * Gets the number of nodes in this graph.
   *
   * @return The number of nodes
   */
  public int getNumNodes() {
    return graph.getNumNodes();
  }

  /**
   * Gets the number of edges in this graph.
   *
   * @return The number of edges.
   */
  public int getNumEdges() {
    return graph.getNumEdges();
  }

  /**
   * Mark a node as being changed.
   *
   * Use this method if the graph weights changed after a previous computation
   * of the max-flow. The next computation will be faster by just considering
   * changed nodes.
   *
   * A node has to be considered changed if any of its adjacent edges changed
   * its weight.
   *
   * @param nodeId The node that changed.
   */
  public void markNode(int node) {

    if (graph.getNextNode(node) == GraphGC.NONE) {
      if (activeQueueLast[1] != GraphGC.NONE)
        graph.setNextNode(activeQueueLast[1], node);
      else
        activeQueueFirst[1] = node;

      activeQueueLast[1] = node;
      graph.setNextNode(node, node);
    }

    graph.isMarked(node, true);
  }

  /*
   * PRIVATE METHODS
   */

  /**
   * Marks a node as being active and adds it to second queue of active nodes.
   */
  private void setNodeActive(int node) {

    if (graph.getNextNode(node) == GraphGC.NONE) {
      if (activeQueueLast[1] != GraphGC.NONE)
        graph.setNextNode(activeQueueLast[1], node);
      else
        activeQueueFirst[1] = node;

      activeQueueLast[1] = node;
      graph.setNextNode(node, node);
    }
  }

  /**
   * Gets the next active node, that is, the first node of the first queue of
   * active nodes. If this queue is empty, the second queue is used. Returns
   * <tt>nyll</tt>, if no active node is left.
   */
  private int getNextActiveNode() {

    int node;

    while (true) {

      node = activeQueueFirst[0];

      if (node == GraphGC.NONE) {
        // queue 0 was empty, try other one
        node = activeQueueFirst[1];

        // swap queues
        activeQueueFirst[0] = activeQueueFirst[1];
        activeQueueLast[0]  = activeQueueLast[1];
        activeQueueFirst[1] = GraphGC.NONE;
        activeQueueLast[1]  = GraphGC.NONE;

        // if other queue was emtpy as well, return Graph.NONE
        if (node == GraphGC.NONE)
          return GraphGC.NONE;
      }

      // remove current node from active list
      if (graph.getNextNode(node) == node) {
        // this was the last one
        activeQueueFirst[0] = GraphGC.NONE;
        activeQueueLast[0]  = GraphGC.NONE;
      } else
        activeQueueFirst[0] = graph.getNextNode(node);

      // not in any list anymore
      graph.setNextNode(node, GraphGC.NONE);

      // return only if it has a parent and is therefore active
      if (graph.getParent(node) != GraphGC.NONE)
        return node;
    }
  }

  /**
   * Mark a node as orphan and add it to the front of the queue.
   */
  private void addOrphanAtFront(int node) {

    graph.setParent(node, GraphGC.ORPHAN);

    orphans.addFirst(node);
  }

  /**
   * Mark a node as orphan and add it to the back of the queue.
   */
  private void addOrphanAtBack(int node) {

    graph.setParent(node, GraphGC.ORPHAN);

    orphans.addLast(node);
  }

  /**
   * Add a node to the list of potentially changed nodes.
   */
  private void addToChangedList(int node) {

    graph.isInChangedList(node, true);
  }

  /**
   * Initialise the algorithm.
   *
   * Only called if <tt>reuseTrees</tt> is false.
   */
  private void maxflowInit() {

    activeQueueFirst[0] = GraphGC.NONE;
    activeQueueLast[0]  = GraphGC.NONE;
    activeQueueFirst[1] = GraphGC.NONE;
    activeQueueLast[1]  = GraphGC.NONE;

    orphans.clear();

    time = 0;

    for (int node = 0; node < graph.getNumNodes(); node++) {

      graph.setNextNode(node, GraphGC.NONE);
      graph.isMarked(node, false);
      graph.isInChangedList(node, false);
      graph.setTimestamp(node, time);

      if (graph.getResidualNodeCapacity(node) > 0) {
        // node is connected to source
        graph.isInSink(node, false);
        graph.setParent(node, GraphGC.TERMINAL);
        setNodeActive(node);
        graph.setDistance(node, 1);
      } else if (graph.getResidualNodeCapacity(node) < 0) {
        // node is connected to sink
        graph.isInSink(node, true);
        graph.setParent(node, GraphGC.TERMINAL);
        setNodeActive(node);
        graph.setDistance(node, 1);
      } else {
        graph.setParent(node, GraphGC.NONE);
      }
    }
  }

  /**
   * Initialise the algorithm.
   *
   * Only called if <tt>reuseTrees</tt> is true.
   */
  private void maxflowReuseTreesInit() {

    int node1;
    int node2;

    int queueStart = activeQueueFirst[1];

    int edge;

    activeQueueFirst[0] = GraphGC.NONE;
    activeQueueLast[0]  = GraphGC.NONE;
    activeQueueFirst[1] = GraphGC.NONE;
    activeQueueLast[1]  = GraphGC.NONE;

    orphans.clear();

    time++;

    while ((node1 = queueStart) != GraphGC.NONE) {

      queueStart = graph.getNextNode(node1);

      if (queueStart == node1)
        queueStart = GraphGC.NONE;

      graph.setNextNode(node1, GraphGC.NONE);
      graph.isMarked(node1, false);
      setNodeActive(node1);

      if (graph.getResidualNodeCapacity(node1) == 0) {
        if (graph.getParent(node1) != GraphGC.NONE)
          addOrphanAtBack(node1);
        continue;
      }

      if (graph.getResidualNodeCapacity(node1) > 0) {

        if (graph.getParent(node1) == GraphGC.NONE || graph.isInSink(node1)) {

          graph.isInSink(node1, false);
          for (edge = graph.getFirstOutgoing(node1); edge != GraphGC.NONE; edge = graph.getNextEdge(edge)) {

            node2 = graph.getHead(edge);
            if (!graph.isMarked(node2)) {
              if (graph.getParent(node2) == graph.getSister(edge))
                addOrphanAtBack(node2);
              if (graph.getParent(node2) != GraphGC.NONE && graph.isInSink(node2) && graph.getResidualEdgeCapacity(edge) > 0)
                setNodeActive(node2);
            }
          }
          addToChangedList(node1);
        }
      } else {

        if (graph.getParent(node1) == GraphGC.NONE || !graph.isInSink(node1)) {

          graph.isInSink(node1, true);
          for (edge = graph.getFirstOutgoing(node1); edge != GraphGC.NONE; edge = graph.getNextEdge(edge)) {

            node2 = graph.getHead(edge);
            if (!graph.isMarked(node2)) {
              if (graph.getParent(node2) == graph.getSister(edge))
                addOrphanAtBack(node2);
              if (graph.getParent(node2) != GraphGC.NONE &&
                  !graph.isInSink(node2) &&
                  graph.getResidualEdgeCapacity(graph.getSister(edge)) > 0)
                setNodeActive(node2);
            }
          }
          addToChangedList(node1);
        }
      }
      graph.setParent(node1, GraphGC.TERMINAL);
      graph.setTimestamp(node1, time);
      graph.setDistance(node1, 1);
    }

    // adoption
    while (orphans.size() > 0) {
      int orphan = orphans.poll();
      if (graph.isInSink(orphan))
        processSinkOrphan(orphan);
      else
        processSourceOrphan(orphan);
    }
  }

  /**
   * Perform the augmentation step of the graph cut algorithm.
   *
   * This is done whenever a path between the source and the sink was found.
   */
  private void augment(int middle) {

    int node;
    int edge;

    float bottleneck;

    // 1. find bottleneck capacity

    // 1a - the source tree
    bottleneck = graph.getResidualEdgeCapacity(middle);
    for (node = graph.getHead(graph.getSister(middle)); ; node = graph.getHead(edge)) {

      edge = graph.getParent(node);

      if (edge == GraphGC.TERMINAL)
        break;
      if (bottleneck > graph.getResidualEdgeCapacity(graph.getSister(edge)))
        bottleneck = graph.getResidualEdgeCapacity(graph.getSister(edge));
    }

    if (bottleneck > graph.getResidualNodeCapacity(node))
      bottleneck = graph.getResidualNodeCapacity(node);

    // 1b - the sink tree
    for (node = graph.getHead(middle); ; node = graph.getHead(edge)) {

      edge = graph.getParent(node);

      if (edge == GraphGC.TERMINAL)
        break;
      if (bottleneck > graph.getResidualEdgeCapacity(edge))
        bottleneck = graph.getResidualEdgeCapacity(edge);
    }
    if (bottleneck > -graph.getResidualNodeCapacity(node))
      bottleneck = -graph.getResidualNodeCapacity(node);

    // 2. augmenting

    // 2a - the source tree
    graph.setResidualEdgeCapacity(graph.getSister(middle), graph.getResidualEdgeCapacity(graph.getSister(middle)) + bottleneck);
    graph.setResidualEdgeCapacity(middle, graph.getResidualEdgeCapacity(middle) - bottleneck);
    for (node = graph.getHead(graph.getSister(middle)); ; node = graph.getHead(edge)) {

      edge = graph.getParent(node);

      if (edge == GraphGC.TERMINAL) {
        // end of path
        break;
      }
      graph.setResidualEdgeCapacity(edge, graph.getResidualEdgeCapacity(edge) + bottleneck);
      graph.setResidualEdgeCapacity(graph.getSister(edge), graph.getResidualEdgeCapacity(graph.getSister(edge)) - bottleneck);
      if (graph.getResidualEdgeCapacity(graph.getSister(edge)) == 0)
        addOrphanAtFront(node);
    }
    graph.setResidualNodeCapacity(node, graph.getResidualNodeCapacity(node) - bottleneck);
    if (graph.getResidualNodeCapacity(node) == 0)
      addOrphanAtFront(node);

    // 2b - the sink tree
    for (node = graph.getHead(middle); ; node = graph.getHead(edge)) {

      edge = graph.getParent(node);

      if (edge == GraphGC.TERMINAL) {
        // end of path
        break;
      }
      graph.setResidualEdgeCapacity(graph.getSister(edge), graph.getResidualEdgeCapacity(graph.getSister(edge)) + bottleneck);
      graph.setResidualEdgeCapacity(edge, graph.getResidualEdgeCapacity(edge) - bottleneck);
      if (graph.getResidualEdgeCapacity(edge) == 0)
        addOrphanAtFront(node);
    }
    graph.setResidualNodeCapacity(node, graph.getResidualNodeCapacity(node) + bottleneck);
    if (graph.getResidualNodeCapacity(node) == 0)
      addOrphanAtFront(node);

    totalFlow += bottleneck;
  }

  /**
   * Adopt an orphan.
   */
  private void processSourceOrphan(int orphan) {

    int bestEdge    = GraphGC.NONE;
    int minDistance = Integer.MAX_VALUE;

    for (int orphanEdge = graph.getFirstOutgoing(orphan); orphanEdge != GraphGC.NONE; orphanEdge = graph.getNextEdge(orphanEdge))
      if (graph.getResidualEdgeCapacity(graph.getSister(orphanEdge)) != 0) {

        int node       = graph.getHead(orphanEdge);
        int parentEdge = graph.getParent(node);

        if (!graph.isInSink(node) && parentEdge != GraphGC.NONE) {

          // check the origin of node
          int distance = 0;
          while (true) {

            if (graph.getTimestamp(node) == time) {
              distance += graph.getDistance(node);
              break;
            }
            parentEdge = graph.getParent(node);
            distance++;
            if (parentEdge == GraphGC.TERMINAL) {
              graph.setTimestamp(node, time);
              graph.setDistance(node, 1);
              break;
            }
            if (parentEdge == GraphGC.ORPHAN) {
              distance = Integer.MAX_VALUE;
              break;
            }
            // otherwise, proceed to the next node
            node = graph.getHead(parentEdge);
          }
          if (distance < Integer.MAX_VALUE) { // node originates from the source

            if (distance < minDistance) {
              bestEdge    = orphanEdge;
              minDistance = distance;
            }
            // set marks along the path
            for (node = graph.getHead(orphanEdge);
                graph.getTimestamp(node) != time;
                node = graph.getHead(graph.getParent(node))) {

              graph.setTimestamp(node, time);
              graph.setDistance(node, distance);
              distance--;
            }
          }
        }
      }

    graph.setParent(orphan, bestEdge);
    if (bestEdge != GraphGC.NONE) {
      graph.setTimestamp(orphan, time);
      graph.setDistance(orphan, minDistance + 1);
    } else {
      // no parent found
      addToChangedList(orphan);

      // process neighbors
      for (int orphanEdge = graph.getFirstOutgoing(orphan); orphanEdge != GraphGC.NONE; orphanEdge = graph.getNextEdge(orphanEdge)) {

        int node = graph.getHead(orphanEdge);
        int parentEdge = graph.getParent(node);
        if (!graph.isInSink(node) && parentEdge != GraphGC.NONE) {

          if (graph.getResidualEdgeCapacity(graph.getSister(orphanEdge)) != 0)
            setNodeActive(node);
          if (parentEdge != GraphGC.TERMINAL && parentEdge != GraphGC.ORPHAN && graph.getHead(parentEdge) == orphan)
            addOrphanAtBack(node);
        }
      }
    }

  }

  /**
   * Adopt an orphan.
   */
  private void processSinkOrphan(int orphan) {

    int bestEdge    = GraphGC.NONE;
    int minDistance = Integer.MAX_VALUE;

    for (int orphanEdge = graph.getFirstOutgoing(orphan); orphanEdge != GraphGC.NONE; orphanEdge = graph.getNextEdge(orphanEdge))
      if (graph.getResidualEdgeCapacity(orphanEdge) != 0) {

        int node       = graph.getHead(orphanEdge);
        int parentEdge = graph.getParent(node);

        if (graph.isInSink(node) && parentEdge != GraphGC.NONE) {

          // check the origin of node
          int distance = 0;
          while (true) {

            if (graph.getTimestamp(node) == time) {
              distance += graph.getDistance(node);
              break;
            }
            parentEdge = graph.getParent(node);
            distance++;
            if (parentEdge == GraphGC.TERMINAL) {
              graph.setTimestamp(node, time);
              graph.setDistance(node, 1);
              break;
            }
            if (parentEdge == GraphGC.ORPHAN) {
              distance = Integer.MAX_VALUE;
              break;
            }
            // otherwise, proceed to the next node
            node = graph.getHead(parentEdge);
          }
          if (distance < Integer.MAX_VALUE) {
            // node originates from the sink
            if (distance < minDistance) {
              bestEdge    = orphanEdge;
              minDistance = distance;
            }
            // set marks along the path
            for (node = graph.getHead(orphanEdge);
                graph.getTimestamp(node) != time;
                node = graph.getHead(graph.getParent(node))) {

              graph.setTimestamp(node, time);
              graph.setDistance(node, distance);
              distance--;
            }
          }
        }
      }

    graph.setParent(orphan, bestEdge);
    if (bestEdge != GraphGC.NONE) {
      graph.setTimestamp(orphan, time);
      graph.setDistance(orphan, minDistance + 1);
    } else {
      // no parent found
      addToChangedList(orphan);

      // process neighbors
      for (int orphanEdge = graph.getFirstOutgoing(orphan); orphanEdge != GraphGC.NONE; orphanEdge = graph.getNextEdge(orphanEdge)) {

        int node = graph.getHead(orphanEdge);
        int parentEdge = graph.getParent(node);
        if (graph.isInSink(node) && parentEdge != GraphGC.NONE) {

          if (graph.getResidualEdgeCapacity(orphanEdge) != 0)
            setNodeActive(node);
          if (parentEdge != GraphGC.TERMINAL && parentEdge != GraphGC.ORPHAN && graph.getHead(parentEdge) == orphan)
            addOrphanAtBack(node);
        }
      }
    }
  }


}
