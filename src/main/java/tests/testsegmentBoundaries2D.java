package tests;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

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
import java.util.concurrent.TimeUnit;

/** An ImageJ2 command with preview capabilities. */ 
@Plugin(type = Command.class,
	menuPath = "Tutorials>Command with Preview")
public class testsegmentBoundaries2D {
	
	
	
	@Parameter 
	static ImagePlus imagePlus;
	
	@Parameter
	static ImageProcessor imageProcessor;

	@Parameter
	static int clusters = 5;
	

	public static void main(String[] args) {	// main-method. starting class.
		final ImageJ ij = net.imagej.Main.launch(args);
		
		int[] labelcolor = {0, 1, 2, 3, 4, 5};
		int[] realcolor = {0, 50, 100, 150, 200, 250};
		int[][] label = generateLabel(labelcolor);
		int [][] image = generateLabel(realcolor);
		
		// Create a beautiful test image.
		long[] dims = {label.length, label[0].length};
		String name = "A test iamge";
		AxisType[] axes = {Axes.X, Axes.Y};
		int bitsPerPixel = 8;
		boolean signed = false;
		boolean floating = false;
		final Dataset dataset =
			ij.dataset().create(dims, name, axes, bitsPerPixel, signed, floating);
		ij.ui().show(dataset);

		imagePlus = IJ.getImage();
		imageProcessor = imagePlus.getProcessor();
		
		for(int x = 0; x < label.length; x++) {
			for(int y = 0; y < label[0].length; y++) {


					imageProcessor.putPixel(x,y, image[x][y]);
			}
		}
		

		ArrayList<ArrayList<int[]>> boundaryCoords2 = segmentBoundaries(label, clusters);
		ArrayList<ArrayList<int[]>> boundariesRaw = segmentBoundariesRaw(label, clusters, CONNECT8); 
		
		//drawborders(boundaryCoords2, imageProcessor, imagePlus);
		drawborders(boundariesRaw, imageProcessor, imagePlus);

	}
	
	
	
	public static int[][] generateLabel(int[] color){
		int[][] label = new int[50][50];
		
		for(int x = 0; x < label.length; x++) {
			for(int y = 0; y < label[0].length; y++) {
					
					if (x < (label.length/2) && y < (label[0].length/2)) {
						label[x][y] = color[0];
					}
					else if (x >= (label.length/2) && y < (label[0].length/2)) {
						label[x][y] = color[1];
					}	
					else if (x < (label.length/2) && y >= (label[0].length/2)) {
						label[x][y] = color[2];
					}					
					else if (x >= (label.length/2) && y >= (label[0].length/2)) {
						label[x][y] = color[3];
					}					

				
			}
		}
		
		for(int x = (label.length/4); x < (3*label.length/4); x++) {
			for(int y = (label[0].length/4); y < (3*label[0].length/4); y++) {
				label[x][y] = color[4];
			}
		}
		
		label[23][8] = color[1];
		label[24][8] = color[1];
		label[20][9] = color[4];
//		label[19][9] = color[4];
//		label[20][8] = color[4];
//		label[19][8] = color[4];
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
	public static ArrayList<ArrayList<int[]>> segmentBoundaries(final int[][] labels, final int nbLabels) {
		// init...
		ArrayList<ArrayList<int[]>> boundaryCoords = new ArrayList<ArrayList<int[]>>();
		
		for (int i=0; i<nbLabels; i++) {
			boundaryCoords.add(i, null );
		}

		// segmentation sizes
		int width = labels.length;
		int height = labels[0].length;
		final int[][] neighbors = CONNECT8;
		boolean bound;
		int label, count, countMax = width*height;
		int k, c, x, y, xT, yT;
		
		// go over all pixels in distance 1 from image boundaries
		for (int i=0; i<width; i++) {
			for (int j=0; j<height; j++) {

				// if this element was not explored
				if (boundaryCoords.get(labels[i][j]) == null) {
					// init the array of coords
					x = i;
					y = j;
					label = labels[x][y];
					count = 0;
					// init boundary segment and add first element
					boundaryCoords.set(label, new ArrayList<int[]>() );
					//boundaryCoords.get(label).add( new int[]{x, y} );
					// direction of last added element
					k = 0;
					// exploring the boundary of the segment by pixel
					do {
						bound = false;
						// over all defined neighbours starting from  previous direction
						for (c=0; c<neighbors.length+1; c++) {
							count ++;
							// get index in bounds of the array
							k = ++k % neighbors.length;
							// temporary coordinates
							xT = x+neighbors[k][0];
							yT = y+neighbors[k][1];
							
							// check if it is inside image
							if (xT<0 || xT>=width || yT<0 || yT>=height) {
								bound = true;
							// if this is not the label set outside the segm.
							} else if (label != labels[ xT ][ yT ]) {
								bound = true;
							// if this is the first point inside the segment
							} else if (bound==true && label==labels[xT][yT]) {
								// add boundary point
								boundaryCoords.get(label).add( new int[]{x, y} );
								x = xT;
								y = yT;
								// next time star in following direction -4
								k += neighbors.length -4;
								break;
							}
						}
					
					// until you come to the first point 
					} while ((x!=boundaryCoords.get(label).get(0)[0] || y!=boundaryCoords.get(label).get(0)[1]) && count<countMax);
										
					// add the initial point again
					boundaryCoords.get(label).add( new int[]{x, y} );
				}
			}
		}
		System.out.println("segmentBoundaries done");
		return boundaryCoords;
	}
	

	/**
	 * returns coordinates of all points belonging to the boundaries among 
	 * different labels in given segmentation, all image boundaries are 
	 * automatically considered also as segments boundaries
	 * Note: there is no order in these boundary points
	 * 
	 * @param neighborhood is matrix int[nbPints][nbDims] describing the 
	 * connectivity around one pixel in given regular grid
	 * @param nbLabels is number of all labelles in segmentation
	 * @return int[nbPixels][nbDims] coordinates of the bordering points
	 */
	public static ArrayList<ArrayList<int[]>> segmentBoundariesRaw(final int[][] labels, final int nbLabels, final int[][] neighbors) {
		ArrayList<ArrayList<int[]>> boundaryCoords = new ArrayList<ArrayList<int[]>>();
		for (int i=0; i<nbLabels; i++) {
			boundaryCoords.add(i, new ArrayList<int[]>() );
		}
		// segmentation sizes
		int width = labels.length;
		int height = labels[0].length;
		int k;

		// add just boundaries because all image boundaries has to be segment boundaries
		// parallel first and last row
		for ( int i = 0; i < width; i++ ) {
			boundaryCoords.get(labels[i][0]).add(new int[]{i,0});
			boundaryCoords.get(labels[i][height-1]).add(new int[]{i,height-1});
		}
		// parallel first and last column
		for ( int i = 1; i < height-1; i++ ) {
			boundaryCoords.get(labels[0][i]).add(new int[]{0,i});
			boundaryCoords.get(labels[width-1][i]).add(new int[]{width-1,i});
		}
		
		// go over all pixels in distance 1 from image boundaries which means 
		// that segm[0][0], segm[0][end], segm[end][0], segm[end][end] may 
		// not be not covered in case of 4-neighborhood
		for (int i=1; i<width-1; i++) {
			for (int j=1; j<height-1; j++) {
				// over all defined neighbors
				for (k=0; k<neighbors.length; k++) {
					if (labels[i][j] != labels[ i+neighbors[k][0] ][ j+neighbors[k][1] ]) {
						boundaryCoords.get(labels[i][j]).add( new int[]{i, j} );
						break;
					}
				}
			}
		}
		
		return boundaryCoords;
	}
	
	
	protected static void drawborders(ArrayList<ArrayList<int[]>> borderCoords, ImageProcessor imageProcessor, ImagePlus imagePlus){

			for (int i = 0; i < clusters; i++) {
				for (int j = 0; j < borderCoords.get(i).size(); j++) {
			
					imageProcessor.putPixel(borderCoords.get(i).get(j)[0], borderCoords.get(i).get(j)[1], 255);
					imagePlus.updateAndDraw();
					
					try {
					    Thread.sleep((1000/imagePlus.getWidth()));
					} catch (InterruptedException e) {
					    // recommended because catching InterruptedException clears interrupt flag
					    Thread.currentThread().interrupt();
					    return;
					}
			}
		}
	}
	
public static final int[][] CONNECT8 = {{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1}};


}


