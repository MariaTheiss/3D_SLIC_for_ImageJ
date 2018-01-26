package tests;
import java.util.ArrayList;
import java.util.Arrays;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;

import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/** An ImageJ2 command with preview capabilities. */ 
@Plugin(type = Command.class,
	menuPath = "Tutorials>Command with Preview")
public class testsegmentBoundaries3D {
	
	
	
	@Parameter 
	static ImagePlus imagePlus;
	
	@Parameter
	static ImageProcessor imageProcessor;
	
	@Parameter
	static ImageStack imageStack;
	
	static int clusters = 8;
	
	public static void main(String[] args) {	// main-method. starting class.
		final ImageJ ij = net.imagej.Main.launch(args);
		
		// labelnumber contains integers that function as label-indices. 
		// realcolor contains the colors of clusters with equivalent indices.
		int[] labelnumber = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
		int[] realcolor = {0, 30, 60, 90, 120, 150, 180, 210, 240};

		// generateLabel creates a 3D array with colors or cluster-indices. 
		// this can be represented as an Image
		int[][][] label = generateLabel(labelnumber);
		int[][][] picture = generateLabel(realcolor);


		// Create a beautiful test image.
		long[] dims = {label.length, label[0].length, label[0][0].length};
		String name = "A test iamge";
		AxisType[] axes = {Axes.X, Axes.Y, Axes.Z};
		int bitsPerPixel = 8;
		boolean signed = false;
		boolean floating = false;
		final Dataset dataset =
			ij.dataset().create(dims, name, axes, bitsPerPixel, signed, floating);
		ij.ui().show(dataset);

		// colorcode the clusters generated with generateLabel 
		imagePlus = IJ.getImage();
		imageStack = imagePlus.getStack();
		
		for(int x = 0; x < picture.length; x++) {
			for(int y = 0; y < picture[0].length; y++) {
				for(int z = 0; z < picture[0][0].length; z++) {

					imageStack.setVoxel(x,y,z, picture[x][y][z]);
				}
			}
		}

		// segmentBoundaries function
	//	ArrayList<ArrayList<int[]>> boundaryCoords = segmentBoundaries2(label, clusters);
		ArrayList<ArrayList<int[]>> boundariesRaw = segmentBoundariesRaw(label, clusters, CONNECT26);
		
		
		// draw Boundaries saved in segmentBoundaries 
	//	drawborders(boundaryCoords, imageProcessor, imagePlus, imageStack);
		drawborders(boundariesRaw, imageProcessor, imagePlus, imageStack);

	//	System.out.println(Arrays.toString(boundaryCoords.toArray()));
		System.out.println("done");
	}
	
	public static int[][][] generateLabel(int[] clustermarking){
		int[][][] label = new int[100][100][100];
		
		for(int x = 0; x < label.length; x++) {
			for(int y = 0; y < label[0].length; y++) {
				for(int z = 0; z < label[0][0].length; z++) {
					
					if (x < (label.length/2) && y < (label[0].length/2) && z < (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[0];
					}
					else if (x >= (label.length/2) && y < (label[0].length/2) && z < (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[1];
					}	
					else if (x < (label.length/2) && y >= (label[0].length/2) && z < (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[2];
					}					
					else if (x >= (label.length/2) && y >= (label[0].length/2) && z < (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[3];
					}					
					else if (x < (label.length/2) && y < (label[0].length/2) && z >= (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[4];
					}					
					else if (x >= (label.length/2) && y < (label[0].length/2) && z >= (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[5];
					}					
					else if (x < (label.length/2) && y >= (label[0].length/2) && z >= (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[6];
					}					
					else if (x >= (label.length/2) && y >= (label[0].length/2) && z >= (label[0][0].length/2)) {
						label[x][y][z] = clustermarking[7];
					}					
					
					
				}
			}
		}
		
		return label;
	}

	
	
	/**
	 * returns coordinates of neighbouring points belonging to the boundaries 
	 * among different labels in given segmentation, it assume that each segment 
	 * can be enclosed by a single close curve and also each segment in the whole 
	 * segmentation has unique label (index)
	 * Note: there is no order in these boundary points
	 * Note2: the chosen connectivity is 8
	 * 
	 * @param neighborhood is matrix int[nbPints][nbDims] describing the 
	 * connectivity around one pixel in given regular grid
	 * @param nbLabels is number of all labelles in segmentation
	 */
	public static ArrayList<ArrayList<int[]>> segmentBoundaries(final int[][][] labels, final int nbLabels) {
		// init...
		ArrayList<ArrayList<int[]>> boundaryCoords = new ArrayList<ArrayList<int[]>>();
		
		for (int i = 0; i < nbLabels; i++) {
			boundaryCoords.add(i, null);
		}

		// segmentation sizes
		int width = labels.length;
		int height = labels[0].length;
		int depth = labels[0][0].length;
//		final int[][] neighbors = CONNECT8;	
		final int[][] neighbors = CONNECT26;	

		boolean bound;
		int label, count, countMax = width * height * depth;	//TODO init count somewhere else- for clarity
		int k, c, x, y, z, xT, yT, zT;
		
		// iterate over whole image 
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int l = 0; l < depth; l++) {
	
					// if this element was not explored
					if (boundaryCoords.get(labels[i][j][l]) == null) {
						// init the array of coords
						
						x = i;
						y = j;
						z = l;

						label = labels[x][y][z];
						
						count = 0;	
						
						// init boundary segment and add first element
						boundaryCoords.set(label, new ArrayList<int[]>() );
						
						// direction of last added element
						k = 0; // Index of neighbor 
						
						// exploring the boundary of the segment by pixel
						do {
							bound = false;
							// over all defined neighbours starting from  previous direction
							for (c = 0; c < neighbors.length+1; c++) {
								count ++;	// to check if iterator is within image boundaries
								// get index in bounds of the array
								
								k = ++k % neighbors.length; 
								
								// temporary coordinates
								xT = x + neighbors[k][0];
								yT = y + neighbors[k][1];
								zT = z + neighbors[k][2];
								
								// check if coordinates are outside of image
								if (xT < 0 || xT >= width || yT < 0 || yT >= height || zT < 0 || zT >= depth) {
									bound = true;
									
								// bound = true if neighbor belongs to a different cluster
								} else if (label != labels[ xT ][ yT ][ zT ]) {
									bound = true;
									
								// If the already visited neighbor[s] belong to a different cluster than 
								// the center-pixel and the currently visited pixel belongs to the same cluster...
								} else if (bound==true && label==labels[xT][yT][zT]) {
									
									// ... add the center-pixel as boundary pixel... 
									boundaryCoords.get(label).add( new int[]{x, y, z} );
									
									// ... and set the currently visited pixel as new center-pixel
									x = xT;
									y = yT;
									z = zT;
									
									// next time star in following direction -4
//									k += neighbors.length -4;	//TODO: change
									k += neighbors.length -13;	//TODO: change

									break;
								}
							}
							// until you come to the first point 
						} while ((	   x != boundaryCoords.get(label).get(0)[0] 
									|| y != boundaryCoords.get(label).get(0)[1] 
									|| z != boundaryCoords.get(label).get(0)[2]) 
									&& count < countMax);
											
						// add the initial point again
						boundaryCoords.get(label).add( new int[]{x, y, z} );
					}
				}
			}
		}
		
		return boundaryCoords;
	}
	

	
	public static ArrayList<ArrayList<int[]>> segmentBoundaries2(final int[][][] labels, final int nbLabels) {
		// init...
		ArrayList<ArrayList<int[]>> boundaryCoords = new ArrayList<ArrayList<int[]>>();
		
		ArrayList<ArrayList<int[]>> boundaryCoords2D = new ArrayList<ArrayList<int[]>>();

		
		for (int i = 0; i < nbLabels; i++) {
			boundaryCoords.add(i, null);
			boundaryCoords2D.add(i, null);
		}

		// segmentation sizes
		int width = labels.length;
		int height = labels[0].length;
		int depth = labels[0][0].length;
		final int[][] neighbors = CONNECT8;	
//		final int[][] neighbors = CONNECT26;	

		boolean bound;
		int label, count, countMax = width * height * depth;	//TODO init count somewhere else- for clarity
		int k, c, x, y, z, xT, yT ;
		
		for (z = 0; z < depth; z++) {
			// iterate over whole image 
			System.out.println("z: " + z);
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
		
		
						// if this element was not explored
						if (boundaryCoords2D.get(labels[i][j][z]) == null) {	//labels[i][j][k] = int (0 - 7) //this is only null for z = 0
							// init the array of coords
							
							x = i;
							y = j;
	
	
							label = labels[x][y][z];
							
							count = 0;	
							
							// init boundary segment and add first element
							if (z < 1){
								boundaryCoords.set(label, new ArrayList<int[]>() );
								boundaryCoords2D.set(label, new ArrayList<int[]>());
							}
							
							// direction of last added element
							k = 0; // Index of neighbor 
							
							// exploring the boundary of the segment by pixel
							do {
								bound = false;
								// over all defined neighbours starting from  previous direction
								for (c = 0; c < neighbors.length+1; c++) {
									count ++;	// to check if iterator is within image boundaries
									// get index in bounds of the array
									
									k = ++k % neighbors.length; 
									
									// temporary coordinates
									xT = x + neighbors[k][0];
									yT = y + neighbors[k][1];
									
									// check if coordinates are outside of image
									if (xT < 0 || xT >= width || yT < 0 || yT >= height || z < 0 || z >= depth) {
										bound = true;
										
									// bound = true if neighbor belongs to a different cluster
									} else if (label != labels[ xT ][ yT ][ z ]) {
										bound = true;
										
									// If the already visited neighbor[s] belong to a different cluster than 
									// the center-pixel and the currently visited pixel belongs to the same cluster...
									} else if (bound==true && label==labels[xT][yT][z]) {
										
										// ... add the center-pixel as boundary pixel... 
										boundaryCoords.get(label).add( new int[]{x, y, z} );
										
										// ... and set the currently visited pixel as new center-pixel
										x = xT;
										y = yT;
										//z = zT;
										
										// next time star in following direction -4
										k += neighbors.length -4;	//TODO: change
										//k += neighbors.length -13;	//TODO: change
	
										break;
									}
								}
								// until you come to the first point 
							} while ((	   x != boundaryCoords.get(label).get(0)[0] 
										|| y != boundaryCoords.get(label).get(0)[1] 
										|| z != boundaryCoords.get(label).get(0)[2]) 
										&& count < countMax);
												
							// add the initial point again
							boundaryCoords.get(label).add( new int[]{x, y, z} );
						}
					}
				}
			}
			
		return boundaryCoords;
	}
	
	
	static void drawborders(ArrayList<ArrayList<int[]>> borderCoords, ImageProcessor imageProcessor, ImagePlus imagePlus, ImageStack imageStack){

			for (int i = 0; i < clusters; i++) {
				for (int j = 0; j < borderCoords.get(i).size(); j++) {
			
					imageStack.setVoxel(borderCoords.get(i).get(j)[0], borderCoords.get(i).get(j)[1], borderCoords.get(i).get(j)[2], 255);
					imagePlus.updateAndDraw();
					
//					try {
//					    Thread.sleep((1000/imagePlus.getWidth()));
//					} catch (InterruptedException e) {
//					    // recommended because catching InterruptedException clears interrupt flag
//					    Thread.currentThread().interrupt();
//					    return;
//					}
			}
		}
	}
	
	public static ArrayList<ArrayList<int[]>> segmentBoundariesRaw(final int[][][] labels, final int nbLabels, final int[][] neighbors) {
		System.out.println("Connectivity3D segmentBoundariesRaw");
		ArrayList<ArrayList<int[]>> boundaryCoords = new ArrayList<ArrayList<int[]>>(); 
		
		for (int i = 0; i < nbLabels; i++) {
			boundaryCoords.add(i, new ArrayList<int[]>() );	//add new array for avery cluster
		}
		// segmentation sizes
		int width = labels.length;
		int height = labels[0].length;
		int depth = labels[0][0].length;	 
		int k;
		
		// add just boundaries because all image boundaries has to be segment boundaries
		//TODO: more efficient loops 
		
		for (int z = 0; z < depth; z++) {
			// parallel first and last row 	
			for ( int i = 0; i < width; i++ ) {
				boundaryCoords.get(labels[i][0][z]).add(new int[]{i, 0, z+1});
				boundaryCoords.get(labels[i][height-1][z]).add(new int[]{i, height-1, z+1});
			}
			// parallel first and last column
			for ( int i = 1; i < height-1; i++ ) {
				boundaryCoords.get(labels[0][i][z]).add(new int[]{0, i, z+1});
				boundaryCoords.get(labels[width-1][i][z]).add(new int[]{width-1, i, z+1});
			}	
		}
		
		// also z 
		for (int i = 1; i < width-1; i++) {
			for (int j = 1; j < height-1; j++) {
				boundaryCoords.get(labels[i][j][0]).add(new int[]{i, j, 1});
				boundaryCoords.get(labels[i][j][depth-1]).add(new int[]{i, j, depth});
			}
		}
		
		
		
		// go over all pixels in distance 1 from image boundaries which means 
		// that segm[0][0], segm[0][end], segm[end][0], segm[end][end] may 
		// not be not covered in case of 4-neighborhood
		for (int i = 1; i < width-1; i++) {
			for (int j = 1; j < height-1; j++) {
				for (int l = 1; l < depth-1; l++) {
					
					// over all defined neighbors
					for (k = 0; k < neighbors.length; k++) {
						
						if (labels[i][j][l] != labels	[ i + neighbors[k][0] ]
														[ j + neighbors[k][1] ]
														[ l + neighbors[k][2] ]) {
							
							boundaryCoords.get(labels[i][j][l]).add( new int[]{i, j, l} );
							
							break;
						}
					}
				}
			}
		}
		
		return boundaryCoords;
	}
	
	
	//TODO: double-check CONNECT26
	public static final int[][] CONNECT26 = {
											{0,	0,	1},	
											{-1, 1,	1},	{-1, 0,	1},	{-1, -1, 1}, {0, -1, 1},	
											{1,	-1,	1},	{1,	0,	1},	{1,	1,	1},	{0,	1,	1},
											{-1, 1, 0}, {-1, 0, 0},{-1 ,-1 ,0},{0 , -1, 0},
											{1, -1, 0}, {1 , 0, 0}, {1, 1, 0}, {0, 1, 0},
											{-1, 1, -1}, {-1, 0, -1},{-1 ,-1 ,-1},{0 , -1, -1},
											{1, -1, -1}, {1 , 0, -1}, {1, 1, -1}, {0, 1, -1},{0, 0, -1}
											};
	
	
	/**
	 * static parameterization of 8-neighbor connectivity 
	 * gives the coordinates on 2D grid
	 */
	public static final int[][] CONNECT8 = {{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1}};
	// public static int[][] CONNECT8 = {{-1,-1,-1,0,1,1,1,0},{1,0,-1,-1,-1,0,1,1}}

}


