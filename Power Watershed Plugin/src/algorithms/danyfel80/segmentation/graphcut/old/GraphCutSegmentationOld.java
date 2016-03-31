/**
 * 
 */
package algorithms.danyfel80.segmentation.graphcut.old;

import icy.image.IcyBufferedImage;
import icy.roi.ROI;
import icy.sequence.Sequence;
import icy.sequence.SequenceUtil;
import icy.type.DataType;
import icy.type.point.Point5D;

import java.awt.Color;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import javax.vecmath.Point3i;

import com.google.common.collect.Lists;

import algorithms.danyfel80.segmentation.SegmentationAlgorithm;


/**
 * @author Daniel Felipe Gonzalez Obando
 * 
 * This class implements the graph cut image segmentation algorithm proposed in
 * "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy 
 * Minimization in Vision" by Yuri Boykov and Vladimir Kolmogorov in IEEE 
 * Transactions on PAMI, Vol. 26, No. 9, pp. 1124-1137, Sept. 2004
 *
 */
public class GraphCutSegmentationOld extends SegmentationAlgorithm {
	
	private Sequence inSequence; // Input sequence
	private Sequence grayInSequence; // Gray-scale input sequence
	private int sizeX; // Size of sequence in X dimension
	private int sizeY; // Size of sequence in Y dimension
	private int sizeZ; // Size of sequence in Z dimension
	//private int sizeC; // Amount of channels of the input sequence

	private Color[] colors; // Colors used for seeds
	private Sequence vertices; // vertices of the graph

	private List<List<Double>> capacities; // edges costs
	private Queue<Integer> active; // active vertices in graph
	private Queue<Integer> orphans; // orphan vertices in graph
	private List<Integer> tree; // indicates to which tree a voxel belongs
	private List<Integer> parent; // parents of the vertices in the graph

	private final double lambda; // lambda constant for weight calibration
	private double K = 0.0; // K is larger than the sum	of all n-links costs for any given pixel p.

	/**
	 * Graph cut constructor. Uses inSequence to execute the segmentation.
	 * @param inSequence Input sequence for the segmentation.
	 */
	public GraphCutSegmentationOld(Sequence inSequence, double lambda) {
		this.lambda = lambda;
		this.inSequence = SequenceUtil.convertToType(inSequence, DataType.DOUBLE, true);
		this.grayInSequence = SequenceUtil.toGray(inSequence);
		this.grayInSequence = SequenceUtil.convertToType(grayInSequence, DataType.DOUBLE, false);

		sizeZ = this.inSequence.getSizeZ();
		sizeX = this.inSequence.getSizeX();
		sizeY = this.inSequence.getSizeY();
		//sizeC = this.inSequence.getSizeC();

		prepareGraph();
	}

	/* (non-Javadoc)
	 * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#prepareGraph()
	 */
	@Override
	public void prepareGraph() {
		vertices = new Sequence(inSequence.getName() + "_segmentation");
		parent = new ArrayList<Integer>((sizeZ*sizeX*sizeY)+2);
		tree = new ArrayList<Integer>((sizeZ*sizeX*sizeY)+2);
		capacities = new ArrayList<List<Double>>((sizeZ*sizeX*sizeY)+2);
		for(int i = 0; i < 2 + (sizeX*sizeY*sizeZ); i++) {
			capacities.add(new ArrayList<Double>((i < 2? sizeX*sizeY*sizeZ: 29))); // capacities for each neighbor node (27 neighbors + S + T)
			tree.add(null);
			parent.add(null);
		}

		double varianceSeq = getIntensityVariance();

		vertices.beginUpdate();
		for (int z = 0; z < sizeZ; z++) {
			vertices.setImage(0, z, new IcyBufferedImage(
					sizeX, sizeY, 1, DataType.DOUBLE));

			//set weights for neighbors of each pixel
			for (int x = 0; x < sizeX; x++) {
				for (int y = 0; y < sizeY; y++) {
					setWeights(x, y, z, varianceSeq);
				}
			}
		}
		vertices.endUpdate();
		System.out.println("K = " + K);
	}

	/**
	 * calculates the intensity variance of inSequence converted to gray scale;
	 * @return intensity variance of input sequence. 
	 */
	private double getIntensityVariance() {
		double[][][] inSequenceData = grayInSequence.getDataXYCZAsDouble(0);
		double meanIp = 0.0;
		double varianceIp = 0.0;
		int amount = 0;

		for (int z = 0; z < sizeZ; z++) {
			for (int x = 0; x < sizeX; x++) {
				for (int y = 0; y < sizeY; y++) {
					meanIp += inSequenceData[z][0][x*sizeY + y];
					amount++;
				}
			}
		}
		meanIp /= (double)amount;

		for (int z = 0; z < sizeZ; z++) {
			for (int x = 0; x < sizeX; x++) {
				for (int y = 0; y < sizeY; y++) {
					double Ip = inSequenceData[z][0][x*sizeY + y] - meanIp;
					varianceIp += Ip*Ip;
				}
			}
		}
		varianceIp /= (double)amount;

		return varianceIp;
	}

	/**
	 * Sets the weights (capacities) for the neighborhood of a specific voxel in the input sequence.
	 * @param x position of voxel in x
	 * @param y position of voxel in y
	 * @param z position of voxel in z
	 * @param varianceSeq intensity variance of input sequence
	 */
	private void setWeights(int x, int y, int z, double varianceSeq) {
		double[][][] inSequenceData = grayInSequence.getDataXYCZAsDouble(0);
		double iP = inSequenceData[z][0][x*sizeY + y]; // Intensity at pixel p.
		List<Double> weightsP = capacities.get(2+ convertXYZToPos(x, y, z));

		weightsP.add(0.0); // dummy weight to S
		weightsP.add(0.0); // dummy weight to T

		double weightsum = 0;

		for (int dz = -1; dz <= 1; dz++) {
			for (int dx = -1; dx <= 1; dx++) {
				for (int dy = -1; dy <= 1; dy++) {
					if (dz != 0 || dx != 0 || dy != 0){
						if (z+dz >= 0 && z+dz < sizeZ && x+dx >= 0 && x+dx < sizeX && y+dy >= 0 && y+dy < sizeY){						
							double iQ = inSequenceData[z+dz][0][(x+dx)*sizeY + (y+dy)];
							double diffI = iP-iQ;
							double weight = Math.exp(-(diffI*diffI)/(2.0*varianceSeq*varianceSeq)) * (1.0/(double)(dx*dx+dy*dy+dz*dz));
							weightsum += weight;
							weightsP.add(weight);
						} else {
							weightsP.add(null);
						}
					} else {
						weightsP.add(null);
					}
				}
			}
		}
		if (weightsum + 1.0 > K) K = weightsum + 1.0;
	}

	/* (non-Javadoc)
	 * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#executeSegmentation(java.util.List)
	 */
	@Override
	public void executeSegmentation(List<ROI> seeds) {
		// Calculate initial segmentation values and capacities.
		List<Color> labels = new ArrayList<>();
		for (ROI roi : seeds) {
			labels.add(roi.getColor());
		}
		colors = labels.toArray(new Color[2]);
		paintInitialSeeds(seeds);
		
		double[][] histograms = getClassesHistograms();
		assignSTCapacities(histograms);

		orphans = new LinkedList<Integer>();
		active = new LinkedList<Integer>();

		tree.set(0, new Integer(0)); // add source source to S
		tree.set(1, new Integer(1)); // add target source to T
		// init active queue
		active.add(0);
		active.add(1);

		// execute max flow
		while(true) {
			List<Integer> path = growSeeds();
			if (path.isEmpty()) {
				break;
			}
			System.out.println(path);
			for (Integer pos : path) {
				if (pos > 1)
					System.out.print(convertPosToXYZ(pos-2) + " -> ");
				else {
					System.out.printf("%d -> ", pos);
				}
			}
			System.out.println("");

			augmentPath(path);
			adoptOrphans();

		}
		System.out.println("end segmentation");

		paintSegmentation();

	}

	/**
	 * Paints the vertices sequence with the seeds color
	 * @param seeds list of seeds (2 rois)
	 */
	private void paintInitialSeeds(List<ROI> seeds) {

		vertices.beginUpdate();
		double[][][] verticesData = vertices.getDataXYCZAsDouble(0);
		for (int z = 0; z < sizeZ; z++) {
			for (int x = 0; x < sizeX; x++) {
				for (int y = 0; y < sizeY; y++) {
					if (seeds.get(0).contains(new Point5D.Double(x,y,z,0,0)))
						verticesData[z][0][x*sizeY + y] = 1;
					else if (seeds.get(1).contains(new Point5D.Double(x,y,z,0,0)))
						verticesData[z][0][x*sizeY + y] = 2;
				}
			}
		}
		vertices.dataChanged();
		vertices.endUpdate();
	}


	/**
	 * Computes the histogram for each of the seed ROIs.
	 * @return Histogram array (one histogram for each ROI).
	 */
	private double[][] getClassesHistograms() {
		int [] amount = new int[2];
		double[][] histograms = new double[2][256];

		double [][][] verticesData = vertices.getDataXYCZAsDouble(0); 
		double [][][] imageData = grayInSequence.getDataXYCZAsDouble(0);

		for (int z = 0; z < sizeZ; z++) {
			for (int x = 0; x < sizeX; x++) {
				for (int y = 0; y < sizeY; y++) {
					double v = verticesData[z][0][x*sizeY + y];
					if (v > 0) {
						int p = ((int)imageData[z][0][x*sizeY + y]);
						histograms[(int)v-1][p] += 1;
						amount[(int)v-1]++;
					}
				}
			}
		}
		for (int i = 0; i < 256; i++) {
			histograms[0][i] /= (double)amount[0];
			histograms[1][i] /= (double)amount[1];
			//System.out.printf("[%d]=%f\n", i, histograms[0][i]);
		}

		return histograms;
	}

	/**
	 * Assigns the capacities from the source and target vertices to each pixel in
	 * the input sequence.
	 * @param histograms histograms of seeds.
	 */
	private void assignSTCapacities(double[][] histograms) {
		double [][][] imageData = grayInSequence.getDataXYCZAsDouble(0);
		double [][][] verticesData = vertices.getDataXYCZAsDouble(0);

		for (int i = 0; i < sizeX*sizeY*sizeZ; i++) {
			 capacities.get(0).add(0.0);
			 capacities.get(1).add(0.0);
		}

		for (int st = 0; st < 2; st++) {
			for (int z = 0; z < sizeZ; z++) {
				for (int x = 0; x < sizeX; x++) {
					for (int y = 0; y < sizeY; y++) {
						int pos = convertXYZToPos(x, y, z);
						if (verticesData[z][0][x*sizeY + y] == st+1) { // if seed
							capacities.get(st).set(pos, K);
							capacities.get((st+1)%2).set(pos, 0.0);
						} else { // if not seed, then use histograms
							int p = (int)imageData[z][0][x*sizeY + y];
							double weight = lambda * (-Math.log(histograms[st][p]));
							if (weight > K) {
								capacities.get(st).set(pos, K);
								capacities.get(pos + 2).set(st, K);
							}
							else {
								capacities.get(st).set(pos, weight);
								capacities.get(pos + 2).set(st, weight);
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Grows the sequence until an augmenting path is found.
	 * @return An augmenting path. Empty list if no path is found.
	 */
	private List<Integer> growSeeds() {
		while (!active.isEmpty()) {

			int posP = active.peek();
			List<Double> neighbors = capacities.get(posP);

			Point3i posPXYZ = null;


			if (posP > 1)
				posPXYZ = convertPosToXYZ(posP-2);

			for (int i = 0; i < neighbors.size(); i++) {
				Double capPQ = neighbors.get(i);
				// Take neighbors with capacity > 0
				if (capPQ != null && capPQ > 0.0) {
					// Set posQ
					int posQ = i;
					if (posP > 1) { // if P is not S or T
						if (posQ > 1) { // if Q is not S or T
							int dz = (i-2)/9 - 1;
							int dx = ((i-2)%9)/3 - 1;
							int dy = (i-2)%3 - 1;
							posQ = convertXYZToPos(posPXYZ.x+dx, posPXYZ.y+dy, posPXYZ.z+dz) + 2;
						} // else Q is either S or T
					} else { // P is either S or T
						posQ += 2;
					}

					// q doesn't belong to a tree
					if (tree.get(posQ) == null) {
						tree.set(posQ, tree.get(posP));
						parent.set(posQ, posP);
						active.add(posQ);
					}
					// if q and p don't belong to same tree -> found and augmenting PATH!!
					if (tree.get(posQ) != null && !tree.get(posP).equals(tree.get(posQ))) {
						List<Integer> path = new ArrayList<Integer>();
						Integer element = posP;
						while (element != null) {
							path.add(element);
							element = parent.get(element);
						}
						path = Lists.reverse(path);
						element = posQ;
						while (element != null) {
							path.add(element);
							element = parent.get(element);
						}
						if (path.get(0) == 1) path = Lists.reverse(path);

						return path;
					}
				}
			}

			active.remove();
		}
		return new ArrayList<Integer>();
	}

	/**
	 * Saturates an augmenting path and finds orphan vertices. 
	 * @param path Augmenting path.
	 */
	private void augmentPath(List<Integer> path) {
		// Find bottleneck capacity
		Double bottleneck = Double.POSITIVE_INFINITY;
		for (int p = 1; p < path.size(); p++) {
			int posP = path.get(p-1);
			int posQ = path.get(p);
			Double capacityPQ = getCapacity(posP, posQ);
//			System.out.println(posP + "->" + posQ + "="+ capacityPQ);
			bottleneck = Math.min(bottleneck, capacityPQ);
		}
//		System.out.println("Bottleneck = " + bottleneck);

		// Push capacity through the path
		if (bottleneck == Double.POSITIVE_INFINITY) { // when all capacities are infinite they exhaust completely
			for (int p = 1; p < path.size(); p++) {
				int posP = path.get(p-1);
				int posQ = path.get(p);
				setCapacity(posP, posQ, 0.0);
				setCapacity(posQ, posP, 0.0);
			}
		}
		else {
			for (int p = 1; p < path.size(); p++) {
				int posP = path.get(p-1);
				int posQ = path.get(p);
				setCapacity(posP, posQ, getCapacity(posP, posQ) - bottleneck);
				setCapacity(posQ, posP, getCapacity(posQ, posP) - bottleneck);
			}
		}

		// Create orphans with exhausted capacities
		for (int p = 1; p < path.size(); p++) {
			int posP = path.get(p-1);
			int posQ = path.get(p);
			Double capacity = getCapacity(posP, posQ);
			if(capacity.equals(0.0)) {
				Integer treeP = tree.get(posP);
				Integer treeQ = tree.get(posQ); 
				if (treeP != null && treeQ != null && treeP.equals(0) && treeQ.equals(0)) { // if P and Q belong tree S
//					System.out.println("0 at "+posP + "->" + posQ);
					parent.set(posQ, null);
					orphans.add(posQ);
				}
				if (treeP != null && treeQ != null && treeP.equals(1) && treeQ.equals(1)) { // if P and Q belong tree T
//					System.out.println("1 at "+posP + "->" + posQ);
					parent.set(posP, null);
					orphans.add(posP);
				}
			}
		}
	}

	/**
	 * Finds parents for orphan vertices in the graph after augmenting a path.
	 */
	private void adoptOrphans() {
		while (!orphans.isEmpty()) {
			// Pick an orphan node
			Integer orp = orphans.poll();

			Point3i orpPos = convertPosToXYZ(orp-2);
//			System.out.println("Orphan:" + orpPos);
			
			// Try to find a valid parent
			boolean nextOrphan = false;
			for (int dz = -1; dz <= 1 && !nextOrphan; dz++) {
				for (int dx = -1; dx <= 1 && !nextOrphan; dx++) {
					for (int dy = -1; dy <= 1 && !nextOrphan; dy++) {
						if(dz != 0 || dx != 0 || dy != 0) {
							if (orpPos.x+dx >= 0 && orpPos.y+dy >= 0 && orpPos.z+dz >= 0 &&
									orpPos.x+dx < sizeX && orpPos.y+dy < sizeY && orpPos.z+dz < sizeZ) {
								Integer neighborPos = convertXYZToPos(orpPos.x+dx, orpPos.y+dy, orpPos.z+dz) + 2;
								Integer origin = neighborPos;
//								System.out.println(origin);
								while(parent.get(origin) != null) {
									origin = parent.get(origin);
								}
//								System.out.println(origin);
								if(tree.get(orp).equals(tree.get(neighborPos)) &&
										getCapacity(neighborPos, orp) > 0.0 &&
										(origin == 0 || origin == 1)) {
									parent.set(orp, neighborPos);
									nextOrphan = true;
								}
							}
						}
					}
				}
			}
			for (int neighborPos = 0; neighborPos < 2 && !nextOrphan; neighborPos++) {
				if(tree.get(orp).equals(tree.get(neighborPos)) &&
						getCapacity(neighborPos, orp) > 0.0) {
					parent.set(orp, neighborPos);
					nextOrphan = true;
				}
			}
			
			// If no parent is found
			if (!nextOrphan) {
				for (int dz = -1; dz <= 1; dz++) {
					for (int dx = -1; dx <= 1; dx++) {
						for (int dy = -1; dy <= 1; dy++) {
							if (dz != 0 || dx != 0 || dy != 0) {
								if (orpPos.x+dx >= 0 && orpPos.y+dy >= 0 && orpPos.z+dz >= 0 &&
										orpPos.x+dx < sizeX && orpPos.y+dy < sizeY && orpPos.z+dz < sizeZ) {
									Integer neighborPos = convertXYZToPos(orpPos.x+dx, orpPos.y+dy, orpPos.z+dz) + 2;
									if (tree.get(orp).equals(neighborPos)) {
										if (getCapacity(neighborPos, orp) > 0.0) {
											active.add(neighborPos);
										}
										if (parent.get(neighborPos).equals(orp)) {
											orphans.add(neighborPos);
											parent.set(neighborPos, null);
										}
									}
								}
							}
						}
					}
				}
				for (int neighborPos = 0; neighborPos < 2 && !nextOrphan; neighborPos++) {
					if (tree.get(orp).equals(tree.get(neighborPos))) {
						if (getCapacity(neighborPos, orp) > 0.0) {
							active.add(neighborPos);
						}
						if (parent.get(neighborPos) != null && parent.get(neighborPos).equals(orp)) {
							orphans.add(neighborPos);
							parent.set(neighborPos, null);
						}
					}
							
				}
				tree.set(orp, null);

				while (active.remove(orp));
			}


		}
	}

	/**
	 * Converts an array position to a XYZ coordinate using input sequence size
	 * @param pos Array position.
	 * @return XYZ coordinate.
	 */
	private Point3i convertPosToXYZ(int pos) {
		Point3i res = new Point3i();
		res.setZ(pos / (sizeY*sizeX));
		res.setX((pos % (sizeY*sizeX)) / sizeY);
		res.setY(pos % sizeY);
		return res;
	}

	/**
	 * Converts a XYZ coordinate to an array position using input sequence size
	 * @param posX X coordinate.
	 * @param posY Y coordinate.
	 * @param posZ Z coordinate.
	 * @return Array position.
	 */
	private int convertXYZToPos(int posX, int posY, int posZ) {
		return posZ*sizeX*sizeY + posX*sizeY + posY;
	}
	
	/**
	 * Gets the capacity between two positions in the graph.
	 * @param posP starting position.
	 * @param posQ target position.
	 * @return capacity P to Q.
	 */
	private Double getCapacity(int posP, int posQ) {
		Double capacity = null;
		if (posP == 0 || posP == 1) {
			capacity = capacities.get(posP).get(posQ-2);
		} else if (posQ == 1 || posQ == 0) {
			capacity = capacities.get(posP).get(posQ);
		} else {
			Point3i pP = convertPosToXYZ(posP-2);
			Point3i pQ = convertPosToXYZ(posQ-2);
			int dx = pQ.x - pP.x;
			int dy = pQ.y - pP.y;
			int dz = pQ.z - pP.z;
			Integer posQrel = (dz+1)*9 + (dx+1)*3 +(dy+1) + 2; 
			capacity = capacities.get(posP).get(posQrel);
		}
		return capacity;
	}
	
	private void setCapacity(int posP, int posQ, double w) {
		if (posP == 0 || posP == 1) {
			capacities.get(posP).set(posQ-2, w);
		} else if (posQ == 1 || posQ == 0) {
			capacities.get(posP).set(posQ, w);
		} else {
			Point3i pP = convertPosToXYZ(posP-2);
			Point3i pQ = convertPosToXYZ(posQ-2);
			int dx = pQ.x - pP.x;
			int dy = pQ.y - pP.y;
			int dz = pQ.z - pP.z;
			Integer posQrel = (dz+1)*9 + (dx+1)*3 +(dy+1) + 2; 
			
			capacities.get(posP).set(posQrel, w);
		}
  }

	private void paintSegmentation() {
		vertices.beginUpdate();
		double[][][] verticesData = vertices.getDataXYCZAsDouble(0);
		for (int z = 0; z < sizeZ; z++) {
			for (int x = 0; x < sizeX; x++) {
				for (int y = 0; y < sizeY; y++) {
					int pos = convertXYZToPos(x, y, z) + 2;
					Integer t = tree.get(pos);
					if (t != null && t == 0) {
						verticesData[z][0][x + y*sizeX] = colors[0].getRed();
						verticesData[z][1][x + y*sizeX] = colors[0].getGreen();
						verticesData[z][2][x + y*sizeX] = colors[0].getBlue();
					}
					else if (t != null && t == 1) {
						verticesData[z][0][x + y*sizeX] = colors[1].getRed();
						verticesData[z][1][x + y*sizeX] = colors[1].getGreen();
						verticesData[z][2][x + y*sizeX] = colors[1].getBlue();
					}
				}
			}
		}
		vertices.dataChanged();
		vertices.endUpdate();
	}

	/* (non-Javadoc)
	 * @see plugins.danyfel80.segmentation.powerwatershed.classes.SegmentationAlgorithm#getSegmentation()
	 */
	@Override
	public Sequence getSegmentation() {
		return vertices;
	}
}
