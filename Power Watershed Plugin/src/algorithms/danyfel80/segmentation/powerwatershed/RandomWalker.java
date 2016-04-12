package algorithms.danyfel80.segmentation.powerwatershed;

/*import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcs;
import edu.emory.mathcs.csparsej.tdouble.Dcs_compress;
import edu.emory.mathcs.csparsej.tdouble.Dcs_lusol;
import edu.emory.mathcs.csparsej.tdouble.Dcs_multiply;
import edu.emory.mathcs.csparsej.tdouble.Dcs_util;*/

import edu.emory.mathcs.csparsej.tfloat.Scs_common.Scs;
import edu.emory.mathcs.csparsej.tfloat.Scs_compress;
import edu.emory.mathcs.csparsej.tfloat.Scs_lusol;
import edu.emory.mathcs.csparsej.tfloat.Scs_multiply;
import edu.emory.mathcs.csparsej.tfloat.Scs_util;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class RandomWalker {
  /**
   * Computes the solution to the Dirichlet problem (RW potential function) 
   * on a general graph represented by an edge list, 
   * given boundary conditions (seeds, etc.).
   * @param edgesIndex List of edges
   * @param M Amount of edges
   * @param verticesIndex List of vertices
   * @param vertexIndicators Boolean array of vertices
   * @param N Number of vertices
   * @param seedsIndex List of nodes that are seeded
   * @param localLabels Associated values for seeds (labels)
   * @param numBoundaries Amount of seeded nodes
   * @param numLabels Amount of different labels
   * @param proba Solution to the Dirichlet problem
   * @return True of random walker is successfully performed, false otherwise.
   */
  public static boolean ExecuteRandomWalker(int[][] edgesIndex, int M,
      int[] verticesIndex, int[] vertexIndicators, int N,
      int[] seedsIndex, float[][] localLabels, int numBoundaries, int numLabels,
      float[][] proba) {
    
    int i, j, k, l, v1, v2;
    boolean[] isSeededVertex = new boolean[N];
    int[] sparseIndicators = new int[N];
    int[] sameEdgeCounts = new int[M];

    // Indexing the edges and the seeds
    for (i = 0; i < N; i++) {
      vertexIndicators[verticesIndex[i]] = i;
    }

    for (j = 0; j < M; j++) {
      v1 = vertexIndicators[edgesIndex[0][j]];
      v2 = vertexIndicators[edgesIndex[1][j]];
      if (v1 < v2) {
        for (i = 0; i < 2; i++) {
          edgesIndex[i][j] = vertexIndicators[edgesIndex[i][j]];
          sparseIndicators[edgesIndex[i][j]]++;
        }
      } else {
        edgesIndex[1][j] = v1;
        edgesIndex[0][j] = v2;
        sparseIndicators[edgesIndex[0][j]]++;
        sparseIndicators[edgesIndex[1][j]]++;
      }
    }
    sortEdges(edgesIndex, M, sameEdgeCounts);

    for (i = 0; i < numBoundaries; i++) {
      seedsIndex[i] = vertexIndicators[seedsIndex[i]];
      isSeededVertex[seedsIndex[i]] = true;
    }

    try {
    // The system to solve is A x = -B X2
    
    Scs A, A2, B, B2;

    // Building matrix A: Laplacian for unseeded nodes
    A2 = Scs_util.cs_spalloc(N-numBoundaries, N-numBoundaries, M*2 + N, true, true);
    if (fillA(A2, N, M, numBoundaries, edgesIndex, isSeededVertex, sparseIndicators, sameEdgeCounts)) {
      // A = compressed-column form of A2
      A = Scs_compress.cs_compress(A2);
      A2 = null;

      // Building boundary matrix B
      B2 = Scs_util.cs_spalloc(N - numBoundaries, numBoundaries, 2*M + N, true, true);
      fillB(B2, N, M, numBoundaries, edgesIndex, isSeededVertex, sparseIndicators, sameEdgeCounts);
      B = Scs_compress.cs_compress(B2);
      B2 = null;

      // Building the right hand side of the system
      Scs X = Scs_util.cs_spalloc(numBoundaries, 1, numBoundaries, true, true);
      Scs X2;
      int rnz, count;
      Scs bTmp;
      float[] b = new float[N - numBoundaries];
      for (l = 0; l < numLabels - 1; l++) {
        // Building vector X
        rnz = 0;
        for (i = 0; i < numBoundaries; i++) {
          X.x[rnz] = localLabels[l][i];
          X.p[rnz] = 0;
          X.i[rnz] = i;
          rnz++;
        }
        X.nz = rnz;
        X.m = numBoundaries;
        X.n = 1;

        X2 = Scs_compress.cs_compress(X);
        bTmp = Scs_multiply.cs_multiply(B, X2);

        for(i = 0; i < N - numBoundaries; i++) {
          b[i] = 0;
        }

        for (i = 0; i < bTmp.nzmax; i++) {
          b[bTmp.i[i]] = -bTmp.x[i];
        }

        // Solve Ax = b by LU decomposition, order = 1
        Scs_lusol.cs_lusol(1, A, b, 1e-7f);

        count = 0;
        for (k = 0; k < N; k++) {
          if (!isSeededVertex[k]) {
            proba[l][verticesIndex[k]] = b[count];
            count++;
          }
        }

        // Enforce boundaries exactly
        for (k = 0; k < numBoundaries; k++) {
          proba[l][verticesIndex[seedsIndex[k]]] = localLabels[l][k];
        }
        X2 = null;
        bTmp = null;
      }

      return true;
    }
    } catch (Exception e) {
      System.err.println("Random walker failed");
      e.printStackTrace();
    }
    return false;
    
  }

  /**
   * Sorts the array of vertices composing edges by ascending node index 
   * and fill an indicator of same edges presence.
   * @param edgesIndex Array of vertices composing edges
   * @param M Amount of edges
   * @param sameEdgeCounts Indicator of same edges presence
   */
  private static void sortEdges(int[][] edgesIndex, int M,
      int[] sameEdgeCounts) {
    int i, j;
    GraphUtils.quickStochasticSortInc(edgesIndex[0], edgesIndex[1], 0, M - 1);
    i = 0;
    while (i < M) {
      j = i;
      while (i < M - 1 && edgesIndex[0][i] == edgesIndex[0][i+1]) {
        i++;
      }
      if (i != j) {
        GraphUtils.quickStochasticSortInc(edgesIndex[1], edgesIndex[0], j, i - j);
      }
      i++;
    }
    
    for (i = 0; i < M; i++) {
      j = 0;
      while (i + j < M - 1 && 
          edgesIndex[0][i + j] == edgesIndex[0][i + j + 1] && 
          edgesIndex[1][i + j] == edgesIndex[1][i + j + 1]) {
        j++;
      }
      sameEdgeCounts[i] = j;
    }
  }

  /**
   * Builds matrix A (Laplacian for unseeded nodes)
   * @param a2 Matrix A to fill
   * @param N Amount of nodes
   * @param M Amount of edges
   * @param numBoundaries Amount of seeds
   * @param edgesIndex Array of node index composing edges
   * @param isSeededVertex Index of seeded nodes
   * @param sparseIndicators Array of index separating seeded and unseeded nodes
   * @param sameEdgeCounts Indicator of same edges presence
   * @return
   */
  private static boolean fillA(Scs a2, int N, int M, int numBoundaries,
      int[][] edgesIndex, boolean[] isSeededVertex, int[] sparseIndicators,
      int[] sameEdgeCounts) {
    int k = 0;
    int rnz = 0;

    // Fill the diagonal
    for (k = 0; k < N; k++) {
      if (!isSeededVertex[k]) {
        a2.x[rnz] = sparseIndicators[k];
        a2.i[rnz] = rnz;
        a2.p[rnz] = rnz;
        rnz++;
      }
    }

    int rnzs = 0;
    int rnzu = 0;

    for (k = 0; k < N; k++) {
      if (isSeededVertex[k]) {
        sparseIndicators[k] = rnzs++;
      } else {
        sparseIndicators[k] = rnzu++;
      }
    }

    for(k = 0; k < M; k++) {
      if (!isSeededVertex[edgesIndex[0][k]] && !isSeededVertex[edgesIndex[1][k]]) {
        a2.x[rnz] = -sameEdgeCounts[k] - 1;
        a2.i[rnz] = sparseIndicators[edgesIndex[0][k]];
        a2.p[rnz] = sparseIndicators[edgesIndex[1][k]];
        rnz++;
        a2.x[rnz] = -sameEdgeCounts[k] - 1;
        a2.p[rnz] = sparseIndicators[edgesIndex[0][k]];
        a2.i[rnz] = sparseIndicators[edgesIndex[1][k]];
        rnz++;
        k += sameEdgeCounts[k];
      }
    }

    a2.nz = rnz;
    a2.m = N - numBoundaries;
    a2.n = N - numBoundaries;
    return true;
  }

  /**
   * Builds matrix B (Laplacian for seeded nodes)
   * @param B Matrix B to fill
   * @param N Amount of nodes
   * @param M Amount of edges
   * @param numBoundaries Amount of seeds
   * @param edgesIndex Array of node index composing edges
   * @param isSeededVertex Index of seeded nodes
   * @param sparseIndicators Array of index separating seeded and unseeded nodes
   * @param sameEdgeCounts Indicator of same edges presence
   */
  private static void fillB(Scs B, int N, int M, int numBoundaries,
      int[][] edgesIndex, boolean[] isSeededVertex, int[] sparseIndicators,
      int[] sameEdgeCounts) {
    int k;
    int rnz;

    rnz = 0;
    for (k = 0; k < M; k++) {
      if (isSeededVertex[edgesIndex[0][k]]) {
        B.x[rnz] = -sameEdgeCounts[k] - 1;
        B.p[rnz] = sparseIndicators[edgesIndex[0][k]];
        B.i[rnz] = sparseIndicators[edgesIndex[1][k]];
        rnz++;
        k += sameEdgeCounts[k];
      } else if (isSeededVertex[edgesIndex[1][k]]) {
        B.x[rnz] = -sameEdgeCounts[k] - 1;
        B.p[rnz] = sparseIndicators[edgesIndex[1][k]];
        B.i[rnz] = sparseIndicators[edgesIndex[0][k]];
        rnz++;
        k += sameEdgeCounts[k];
      } 
    }

    B.nz = rnz;
    B.m = N - numBoundaries;
    B.n = numBoundaries;
  }
}
