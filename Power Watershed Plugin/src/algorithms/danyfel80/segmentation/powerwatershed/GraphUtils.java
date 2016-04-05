package algorithms.danyfel80.segmentation.powerwatershed;

import java.util.ArrayDeque;
import java.util.Deque;

import icy.util.Random;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class GraphUtils {
  /**
   * Compute the index of the j_th edge neighbor of the node "i". 
   *   4
   * 3 0 1   5 front slice,  6 back slice
   *   2 
   * return -1 if the neighbor is outside the image
   * @param i Index of source node
   * @param j Number of the desired edge neighbor of node i
   * @param sx Size of sequence in X
   * @param sy Size of sequence in Y
   * @param sz Size of sequence in Z
   * @return Index of the j_th edge neighbor of the node "i". -1 if the neighbor is outside the image
   */
  public static int getNeighborNodeEdge(int i, int j, int sx, int sy, int sz) {
    int sxy = sx*sy;
    int z = i/sxy;
    int zp = i%sxy;

    int V = (sy-1)*sx;
    int H = (sx-1)*sy;

    switch (j) {
    case 1:
      if (zp%sx >= sx - 1) {
        return -1;
      } else {
        return (zp + V) - (zp/sx) + z*(V + H + sxy);
      }
    case 2:
      if (zp/sx >= sy - 1) {
        return -1;
      } else {
        return zp + z*(V + H + sxy);
      }
    case 3:
      if (zp%sx == 0) {
        return -1;
      } else {
        return (zp + V) - (zp/sx) - 1 + z*(V + H + sxy);
      }
    case 4:
      if (zp/sx == 0) {
        return -1;
      } else {
        return zp - sx + z*(V + H + sxy);
      }
    case 5:
      if (z == 0) {
        return -1;
      } else {
        return z*(V + H) + zp + (z - 1) * sxy;
      }
    case 6:
      if (z >= sz - 1) {
        return -1;
      } else {
        return (z + 1)*(V + H) + zp + z*sxy;
      }
    case -1:
      return i + sxy;
    case 0:
      return i + sxy*2;
    }
    return -1; // should never happen
  }

  /**
   * Computes the index of the k_th neighbor 
   * (only works in 2D, a little faster than the neighbor_edge3D)
   *       1       _|_          
   *     2 0 6     _|_          2 1      _|_|_
   *     3 0 5      |         3 0 0 6     | | 
   *       4                    4 5 
   *                            
   * return -1 if the neighbor is outside the image
   * 
   * indexing edges 2D
   * 
   *     _4_ _5_
   *   0|  1|
   *    |_6_|_7_
   *   2|  3|
   *    |   |
   * 
   * @param i Edge index.
   * @param k Number of the desired neighbor of edge i. 
   * @param sx Size of sequence in X.
   * @param sy Size of sequence in Y.
   * @param sz Size of sequence in Z.
   * @return The index of the k_th neighbor 
   */
  public static int getNeighborEdge(int i, int k, int sx, int sy, int sz) {
    int V = (sy-1)*sx; // nb vertical edges
    
    if (i >= V) { //horizontal
      switch(k)
      {
      case 2:
        if ((i-V) < sx-1) {
          return -1;
        } else {
          return ((i-V)/(sx-1)-1)*sx + ((i - V)%(sx-1));
        }
      case 3:
        if ((i-V)%(sx-1)==0) {
          return -1;
        } else {
          return i-1;
        }
      case 4:
        if (i>(sx-1)*sy+ V -sx){
          return -1;
        } else {
          return ((i-V)/(sx-1)-1)*sx + ((i - V)%(sx-1)) + sx;
        }
      case 5:
        if (i>(sx-1)*sy+ V - sx) {
          return -1;
        } else return ((i-V)/(sx-1)-1)*sx + ((i - V)%(sx-1))  + sx +1;
      case 6:
        if ((i-V)%(sx-1)==sx-2){
          return -1;
        } else {
          return i+1;
        }
      case 1:
        if (i-V<sx-1) {
          return -1;
        } else {
          return ((i-V)/(sx-1)-1)*sx + ((i - V)%(sx-1))+1;
        }
      }
    }
    else { //vertical
      switch(k)
      {
      case 6:
        if (i %sx == sx-1){
          return -1;
        } else {
          return (i+V)-(i/sx);
        }
      case 1:
        if (i < sx) {
          return -1;
        } else {
          return i-sx;
        }
      case 2:
        if (i%sx==0) {
          return -1;
        } else {
          return (i+V)-(i/sx)-1;
        }
      case 3:
        if (i%sx==0) {
          return -1;
        } else {
          return (i+V)-(i/sx)-1+sx-1;
        }
      case 4:
        if (i>=V-sx) {
          return -1;
        } else {
          return i+sx;
        }
      case 5:
        if (i %sx == sx-1) {
          return -1;
        } else {
          return (i+V)-(i/sx)+sx-1;
        }
      }
    }
    return -1; //never happens 
  }
  
  /**
   * Computes the index of the k_th neighbor
   *  
   *      3       _|_       
   *  1 2 0 4 5   _|_       
   *  6 7 0 9 10   |        
   *      8                 
   *                        
   * return -1 if the neighbor is outside the image
   *
   *
   * example of indexing of the edges in 3D
   *
   *        27      28
   *       _____ ________
   *     /      /       /|
   *  12/    13/     14/ |
   *   /__6___/___7__ /  | 23
   *   |      |      |   |
   *  0|     1|     2|17/|
   *   |      |      | / |
   *   |__8___|___9__|/  | 26
   *   |      |      |   |
   *  3|     4|     5|20/
   *   |      |      | /
   *   |__10__|__11__|/
   *   
   * @param node1 Index of node 1.
   * @param node2 Index of node 2.
   * @param i Edge index.
   * @param k Number of the desired neighbor of edge i.
   * @param sx Size of sequence in X.
   * @param sy Size of sequence in Y.
   * @param sz Size of sequence in Z.
   * @return The index of the k_th neighbor.
   */
  public static int getNeighborEdge3D(int node1, int node2, int i, int k, int sx, int sy, int sz) {
    if (sz == 1) {
      return getNeighborEdge(i, k, sx, sy, sz); 
    }
    
    int index= -1;
    int sxy = sx*sy;
    int V = (sy - 1)*sx;
    int H = (sx - 1)*sy;
    if(k <= 6) {
      int zp = node1%(sxy); 
      int z = node1/(sxy);   
      switch(k) {
      case 1:
        if (zp%sx >= sx - 1) {
          return -1;
        } else {
          index = (zp + V) - (zp/sx) + z*(V + H + sxy);
        }
        break;
      case 3:
        if (zp%sx == 0) {
          return -1;
        } else {
          index = (zp + V) - (zp/sx) - 1 + z*(V + H + sxy);
        }
        break;
      case 2:
        if (zp/sx >= sy - 1) {
          return -1;
        }
        else {
          index = zp + z*(V + H + sxy);
        }
        break;
      case 4:
        if (zp/sx == 0) {
          return -1;
        }
        else {
          index = zp - sx + z*(V + H + sxy);
        }
        break;
      case 5:
        if (z == 0) {
          return -1;
        } else {
          index = z*(V + H) + zp + (z - 1)*sxy;
        }
        break;
      case 6:
        if (z >= sz - 1) {
          return -1;
        } else {
          index = (z + 1)*(V + H) + zp + z*sxy;
        }
        break;
      }
    }
    else {
      int zp = node2%(sxy); 
      int z = node2/(sxy);   
      switch(k - 6) {
      case 1:
        if (zp%sx >= sx - 1) {
          return -1;
        } else {
          index = (zp + V) - (zp / sx) + z*(V + H + sxy);
        }
        break;
      case 3:
        if (zp%sx == 0) {
          return -1;
        }
        else {
          index = (zp + V) - (zp/sx) - 1 + z*(V + H + sxy);
        }
        break;
      case 2:
        if (zp/sx >= sy - 1){
          return -1;
        }
        else {
          index = zp + z*(V + H + sxy);
        }
        break;
      case 4:
        if (zp/sx == 0) {
          return -1;
        }
        else {
          index = zp - sx + z*(V + H + sxy);
        }
        break;
      case 5:
        if (z == 0) {
          return -1;
        }
        else {
          index = z*(V + H) + zp + (z - 1)*sxy;
        }
        break;
      case 6:
        if (z >= sz - 1) {
          return -1;
        }
        else {
          index = (z + 1)*(V + H) + zp + z*sxy;
        }
        break;
      }
    }
    
    if (index == i) {
      return -1;
    }
    return index;
  }
  
  
  /**
   * Sort the seedFnSize values of array seedsFunction by ascending order and stocks result in sortedEdges
   * @param seedsFunction Array to sort
   * @param sortedEdges Array such as I[i]=i
   * @param seedsFnSize Number of elements to sort
   * @param numBins Number of buckets
   */
  public static void sort(int[] seedsFunction, int[] sortedEdges, int seedsFnSize, int numBins) {
    int i, j;
    int[] H = new int[numBins];
    for(i = 0; i < seedsFnSize; i++) {
      H[seedsFunction[i]]++;
    }
    @SuppressWarnings("unchecked")
    Deque<Integer>[] buckets =  new Deque[numBins];
    for (i = 0; i < numBins; i++) {
      buckets[i] = new ArrayDeque<Integer>(H[i]);
    }

    for (i = 0; i < seedsFnSize; i++) {
      buckets[seedsFunction[i]].addFirst(i);
    }

    j = 0;
    for (i = numBins - 1; i >= 0; i--) {
      while(!buckets[i].isEmpty()) {
        sortedEdges[j] = buckets[i].removeFirst();
        seedsFunction[j] = i;
        j++;
      }
    }
  }
  
  /**
   * Sorts array values in A from index p (inclusive) to index r (inclusive)
   * in ascending order. 
   * @param A
   * @param I
   * @param p
   * @param r
   */
  public static void quickStochasticSort(int[] A, int[] I, int p, int r) {
    int q;
    if (p < r) {
      q = stochasticPartition(A, I, p, r);
      quickStochasticSort(A, I, p, q);
      quickStochasticSort(A, I, q+1, r);
    }
  }
  
  /**
   * Finds the partition index of A between index p (inclusive) and 
   * index r (inclusive). First partition where elements <= A[q] and second
   * with the other elements (q is randomly selected in the range [p, r]).
   * @param A
   * @param I
   * @param p
   * @param r
   * @return
   */
  private static int stochasticPartition(int[] A, int[] I, int p, int r) {
    int t, t1, q;
    q = p + (Random.nextInt(r - p + 1));
    t = A[p];
    A[p] = A[q];
    A[q] = t;
    
    t1 = I[p];
    I[p] = I[q];
    I[q] = t1;
    
    return findPartition(A, I, p, r);
  }

  /**
   * Finds the partition index of the elements in A between p (inclusive) 
   * and r (inclusive) to split in two groups: those <= A[p] and other 
   * elements.
   *  
   * @param A
   * @param I
   * @param p
   * @param r
   * @return
   */
  private static int findPartition(int[] A, int[] I, int p, int r) {
    
    int  t, t1;
    int x = A[p];
    int i = p - 1;
    int j = r + 1;
    
    while (true) {
      do {
        j--;
      } while (A[j] > x);
      
      do {
        i++;
      } while (A[i] < x);
      
      if (i < j) { 
        t = A[i];
        A[i] = A[j];
        A[j] = t; 
        t1 = I[i];
        I[i] = I[j];
        I[j] = t1; 
      }
      else {
        return j;
      }
    }   
  }

  /**
   * Finds the farthest ancestor node of an specified node. 
   * @param childEdge
   * @param fathers
   * @return the root parent node of the specified child node.
   */
  public static int findElement(int childEdge, int[] fathers) {
    if (fathers[childEdge] != childEdge) {
      fathers[childEdge] = findElement(fathers[childEdge], fathers);
    }
    return fathers[childEdge];
  }

  /**
   * Link two edges as parent and child edges according to their current ranks.
   * @param x Edge 1
   * @param y Edge 2
   * @param ranks Edges' rank array
   * @param fathers Edges' father index
   * @return The father edge between x and y
   */
  public static int linkElements(int x, int y, int[] ranks,
      int[] fathers) {
    if (ranks[x] > ranks[y]) {
      int tmp = x;
      x = y;
      y = tmp;
    }
    if (ranks[x] == ranks[y]) {
      ranks[y] = ranks[y] + 1;
    }
    fathers[x] = y;
    return y;
  }

  
}
