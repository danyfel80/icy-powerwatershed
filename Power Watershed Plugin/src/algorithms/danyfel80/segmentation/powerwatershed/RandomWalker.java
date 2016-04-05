package algorithms.danyfel80.segmentation.powerwatershed;

import edu.emory.mathcs.csparsej.tdouble.Dcs_common.Dcs;
import edu.emory.mathcs.csparsej.tdouble.Dcs_compress;
import edu.emory.mathcs.csparsej.tdouble.Dcs_lusol;
import edu.emory.mathcs.csparsej.tdouble.Dcs_multiply;
import edu.emory.mathcs.csparsej.tdouble.Dcs_util;

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
   * @param vertexIndex List of vertices
   * @param verticesIndicators Boolean array of vertices
   * @param N Number of vertices
   * @param seedsIndex List of nodes that are seeded
   * @param boundaryValues Associated values for seeds (labels)
   * @param numBoundaries Amount of seeded nodes
   * @param numLabels Amount of different labels
   * @param proba Solution to the Dirichlet problem
   * @return True of random walker is successfully performed, false otherwise.
   */
  public static boolean ExecuteRandomWalker(int[][] edgesIndex, int M,
      int[] vertexIndex, int[] verticesIndicators, int N,
      int[] seedsIndex, double[][] boundaryValues, int numBoundaries, int numLabels,
      double[][] proba) {

    int i, j, k, l, v1, v2;
    boolean[] seededVertices = new boolean[N];
    int[] sparseIndicators = new int[N];
    int[] sameEdgeCounts = new int[M];

    // Indexing the edges and the seeds
    for (i = 0; i < N; i++) {
      verticesIndicators[vertexIndex[i]] = i;
    }

    for (j = 0; j < M; j++) {
      v1 = verticesIndicators[edgesIndex[0][j]];
      v2 = verticesIndicators[edgesIndex[1][j]];
      if (v1 < v2) {
        for (i = 0; i < 2; i++) {
          edgesIndex[i][j] = vertexIndex[edgesIndex[j][i]];
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
      seedsIndex[i] = verticesIndicators[seedsIndex[i]];
      seededVertices[seedsIndex[i]] = true;
    }

    // The system to solve is A x = -B X2
    Dcs A, A2, B, B2;

    // Building matrix A: Laplacian for unseeded nodes
    A2 = Dcs_util.cs_spalloc(N-numBoundaries, N-numBoundaries, M*2 + N, true, true);
    if (fillA(A2, N, M, numBoundaries, edgesIndex, seededVertices, sparseIndicators, sameEdgeCounts)) {
      // A = compressed-column form of A2
      A = Dcs_compress.cs_compress(A2);

      // Building boundary matrix B
      B2 = Dcs_util.cs_spalloc(N - numBoundaries, numBoundaries, 2*M + N, true, true);
      fillB(B2, N, M, numBoundaries, edgesIndex, seededVertices, sparseIndicators, sameEdgeCounts);
      B = Dcs_compress.cs_compress(B2);

      // Building the right hand side of the system
      Dcs X = Dcs_util.cs_spalloc(numBoundaries, 1, numBoundaries, true, true);
      Dcs X2;
      int rnz, count;
      Dcs bTmp;
      double[] b = new double[N - numBoundaries];
      for (l = 0; l < numLabels - 1; l++) {
        // Building vector X
        rnz = 0;
        for (i = 0; i < numBoundaries; i++) {
          X.x[rnz] = boundaryValues[l][i];
          X.p[rnz] = 0;
          X.i[rnz] = i;
          rnz++;
        }
        X.nz = rnz;
        X.m = numBoundaries;
        X.n = 1;

        X2 = Dcs_compress.cs_compress(X);
        bTmp = Dcs_multiply.cs_multiply(B, X2);

        for(i = 0; i < N - numBoundaries; i++) {
          b[i] = 0;
        }

        for (i = 0; i < bTmp.nzmax; i++) {
          b[bTmp.i[i]] = -bTmp.x[i];
        }

        // Solve Ax = b by LU decomposition, order = 1
        Dcs_lusol.cs_lusol(1, A, b, 1e-7);

        count = 0;
        for (k = 0; k < N; k++) {
          if (seededVertices[k] == false) {
            proba[l][vertexIndex[k]] = b[count];
            count++;
          }
        }

        // Enforce boundaries exactly
        for (k = 0; k < numBoundaries; k++) {
          proba[l][vertexIndex[seedsIndex[k]]] = boundaryValues[l][k];
        }
      }

      return true;
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
    GraphUtils.quickStochasticSort(edgesIndex[0], edgesIndex[1], 0, M - 1);
    i = 0;
    while (i < M) {
      j = i;
      while (i < M - 1 && edgesIndex[0][i] == edgesIndex[0][i+1]) {
        i++;
      }
      if (i != j) {
        GraphUtils.quickStochasticSort(edgesIndex[1], edgesIndex[0], j, i - j);
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
   * @param A Matrix A to fill
   * @param N Amount of nodes
   * @param M Amount of edges
   * @param numBoundaries Amount of seeds
   * @param edgesIndex Array of node index composing edges
   * @param seededVertices Index of seeded nodes
   * @param sparseIndicators Array of index separating seeded and unseeded nodes
   * @param sameEdgeCounts Indicator of same edges presence
   * @return
   */
  private static boolean fillA(Dcs A, int N, int M, int numBoundaries,
      int[][] edgesIndex, boolean[] seededVertices, int[] sparseIndicators,
      int[] sameEdgeCounts) {
    int k = 0;
    int rnz = 0;

    // Fill the diagonal
    for (k = 0; k < N; k++) {
      if (seededVertices[k] == false) {
        A.x[rnz] = sparseIndicators[k];
        A.i[rnz] = rnz;
        A.p[rnz] = rnz;
        rnz++;
      }
    }

    int rnzs = 0;
    int rnzu = 0;

    for (k = 0; k < N; k++) {
      if (seededVertices[k]) {
        sparseIndicators[k] = rnzs;
        rnzs++;
      } else {
        sparseIndicators[k] = rnzu;
        rnzu++;
      }
    }

    for(k = 0; k < M; k++) {
      if (!seededVertices[edgesIndex[0][k]] && !seededVertices[edgesIndex[1][k]]) {
        A.x[rnz] = -sameEdgeCounts[k] - 1;
        A.i[rnz] = sparseIndicators[edgesIndex[0][k]];
        A.p[rnz] = sparseIndicators[edgesIndex[1][k]];
        rnz++;
        A.x[rnz] = -sameEdgeCounts[k] - 1;
        A.i[rnz] = sparseIndicators[edgesIndex[0][k]];
        A.p[rnz] = sparseIndicators[edgesIndex[1][k]];
        rnz++;
        k = k + sameEdgeCounts[k];
      }
    }

    A.nz = rnz;
    A.m = N - numBoundaries;
    A.n = N - numBoundaries;
    return true;
  }

  /**
   * Builds matrix B (Laplacian for seeded nodes)
   * @param B Matrix B to fill
   * @param N Amount of nodes
   * @param M Amount of edges
   * @param numBoundaries Amount of seeds
   * @param edgesIndex Array of node index composing edges
   * @param seededVertices Index of seeded nodes
   * @param sparseIndicators Array of index separating seeded and unseeded nodes
   * @param sameEdgeCounts Indicator of same edges presence
   */
  private static void fillB(Dcs B, int N, int M, int numBoundaries,
      int[][] edgesIndex, boolean[] seededVertices, int[] sparseIndicators,
      int[] sameEdgeCounts) {
    int k;
    int rnz;

    rnz = 0;
    for (k = 0; k < M; k++) {
      if (seededVertices[edgesIndex[0][k]]) {
        B.x[rnz] = -sameEdgeCounts[k] - 1;
        B.p[rnz] = sparseIndicators[edgesIndex[0][k]];
        B.i[rnz] = sparseIndicators[edgesIndex[1][k]];
        rnz++;
        k = k + sameEdgeCounts[k];
      } else if (seededVertices[edgesIndex[1][k]]) {
        B.x[rnz] = -sameEdgeCounts[k] - 1;
        B.p[rnz] = sparseIndicators[edgesIndex[1][k]];
        B.i[rnz] = sparseIndicators[edgesIndex[0][k]];
        rnz++;
        k = k + sameEdgeCounts[k];
      } 
    }

    B.nz = rnz;
    B.m = N - numBoundaries;
    B.n = numBoundaries;
  }
}
