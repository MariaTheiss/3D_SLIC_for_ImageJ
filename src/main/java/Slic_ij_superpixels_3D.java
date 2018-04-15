package jSLIC3D;
/**
 * @file
 */

import java.awt.Color;

import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.StackProcessor;

import java.lang.Math;
import jSLIC3D.ThreadAssignment;
import jSLIC3D.ConvertImage;
import jSLIC3D.ThreadUpdate;
import jSLIC3D.Threading;
import jSLIC3D.Labelling3D;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import imagescience.feature.Edges;
import imagescience.image.Image;

import net.imagej.ImageJ;


@Plugin(type = Command.class,
	menuPath = "Plugins>SLIC_superpixels_3D")
public class Slic_ij_superpixels_3D  implements Command {
	
	

	// segmentation paramters
	protected int gridSize;


	// -- Variables --

	private int widthx, heighty, depthz;		// Dimensions of image stack [px]. Only for 3D images
	private int[][][][] img3D = null;	// original image converted into LAB color space dim int[Width][Height][Depth][Channels]
	
	private int[][] clusterPosition = null;		// vector of cluster-center positions - dim int[nbClusters][n dims (3 for 3D)]
	private int[][] clusterColor = null;		// vector of clusters color - dim int[nbClusters][channels]. For 8 bit gv channels = 1
												// ... clusterColor is the mean L, A and B value for each cluster. 
	
	private int[][][] labels3D = null;			// labeling for each image pixel - dim int[Widthx][Heighty][Depthz]. Number corresponds to cluster.
	private float[][][] distances3D = null;		// minimal distance per pixel according to assigned label (= cluster) 
	private int nbLabels; 						// number of estimated segments (labels)

	private float errThreshold = 0.1f;			// stopping threshold in percent of initial error
	private float factor;						// factor results from @parameter regul. Used to weight distance vs. color measure.
	private float[] distGrid = null;			// distances of each pixel in a ~2 * 2 * 2 gridsize searchwindow from the cluster center. Array is 1D (z -> y -> x)

	Color clrOverlay = null; 

	
	// set of colours
//	private Map<String, Color> colorMap= new HashMap<String, Color>();
	private Map<String, Color> colorMap= new HashMap<String, Color>();
	private Color[] clrs = {null, Color.WHITE, Color.YELLOW, Color.RED, Color.GREEN, Color.BLUE, Color.BLACK};	//FIXME
	private String[] clrString  = {"none", "white", "yellow", "red", "green", "blue", "black"};

	// -- Parameter --
	@Parameter
	ImagePlus imagePlus;	

	@Parameter
	ImageStack stack;
	
	@Parameter(label = "Grid size [px]")
	int gSize = 40;
	
	@Parameter(min = "0f", stepSize = "0.1f", label = "Regul.")	//TODO: check if this is just in the GUI, or the real step size
	float regul = 0.5f;
	
	@Parameter(label = "Show ROIs")
	boolean showROIs = true;

	@Parameter (choices  = {"none", "white", "yellow", "red", "green", "blue", "black"}, label = "Border overlay color")
	String color;
	

	public void initOverlayColor() {
		
		for (int i = 0; i < clrs.length; i++) {
			colorMap.put(clrString[i], clrs[i]);
		}
		System.out.println(colorMap.get(color));
		clrOverlay = colorMap.get(color);
		
	}
	

	// --Command methods --
	@Override
	public void run() {
		
		long startTime, estimTime;
		startTime = System.currentTimeMillis();
		
		try {
			initOverlayColor();
			imagePlus.lock();	
			process(gSize, regul, 9, 0.1f);	
			imagePlus.unlock();
		} catch (Exception e) {
			System.out.println("jSLIC_superpixels_3D2 run: an exception has been thrown: " + e);
		} finally {
			System.out.println("jSLIC_superpixels_3D2 run unlock image");
			imagePlus.unlock();
			estimTime = System.currentTimeMillis() - startTime;

			printInfo("runtime is: " + estimTime/1000 + " s");
		}
		
	}

	
	/**
	 * the whole processing method
	 */
	public void process (int grid, float reg, int maxIter, float sizeThreshold) {
		long startTime, estimTime;
		startTime = System.currentTimeMillis();
		
		IJ.showProgress(0.);
		//System.out.println("initVariables");
		initVaribales(grid, reg);	// init jSLIC superpixels			
		convertImage();
	
		IJ.showProgress(20.);

		printInfo("Processing. Running with gridsize: " + Integer.toString(gridSize) + ", regularity " + Float.toString(regul) + ", factor: " + ((regul * regul) * (float)(gridSize)));
		
		initClusters();
		
		float err, lastErr = Float.MAX_VALUE;
//		float initErr = computeResidualError();	
//		System.out.println(initErr);
		float initErr = 0f;
	
		computeDistGrid();

		for (int i = 0; i < maxIter; i++) {
			//for (int i = 0; i < 2; i++) {
		
			//IJ.showProgress(((i+1)/(double)maxIter) * 100.);
			
			assignmentParallel();

			err = computeResidualError();
			printInfo("Iter " + Integer.toString(i + 1) + ", inter. distance is " + Float.toString(err));
			
			updateParallel();
			System.out.println(Float.MAX_VALUE + 10000000000000f);

			// STOP if consecutive errors are smaller then given threshold
			if ( (lastErr - err) < (initErr * errThreshold)) {	//initErr * errThreshold is always 0
				printInfo("Terminated with diff error " + (lastErr-err));
				i = maxIter;
			} else {
				lastErr = err;
			}
			
		}

		
		int maxlabel = 0;
		for (int i = 0; i < labels3D.length; i++) {
			for (int j = 0; j < labels3D[0].length; j++) {
				for (int k = 0; k < labels3D[0][0].length; k++) {
					
					if(maxlabel < labels3D[i][j][k]) {
						maxlabel = labels3D[i][j][k];
					}
				}
			}
		}
		
		printInfo("number of labels: " + maxlabel);
	

		enforceLabelConnectivity();
//		printInfo("number of labels: " + nbLabels);

		
//		maxlabel = 0;
//		for (int i = 0; i < labels3D.length; i++) {
//			for (int j = 0; j < labels3D[0].length; j++) {
//				for (int k = 0; k < labels3D[0][0].length; k++) {
//					
//					if(maxlabel < labels3D[i][j][k]) {
//						maxlabel = labels3D[i][j][k];
//					}
//				}
//			}
//		}
//		
//		printInfo("number of labels: " + (maxlabel + 1));

		
		
		IJ.showProgress(90.);
		estimTime = System.currentTimeMillis() - startTime;	
		printInfo("processing took: " + estimTime/1000 + " s");
			
		showSegmentation ();
					
		printInfo("jSLIC3D finished.");
		IJ.showProgress(100.);
		


	}
	
	
	protected void initVaribales(int grid, float reg) {
		printInfo("jSLIC3D initialisation...");
		
		long startTime, estimTime;
		startTime = System.currentTimeMillis();

		this.imagePlus = IJ.getImage();	
		this.stack = imagePlus.getStack();
		
		// Check dimensions
		int[] dims = imagePlus.getDimensions();	// dim int[width][height][nChannels][nSlices][nFrames]
		if ((dims[2] != 1) || (dims[3] * dims[4]) == 1){	// TODO: read in frames as depthz if depth = 1 and frames != 1
			IJ.error("ERROR: Unsupported image dimensions!");
		}
		
		// Get dimensions
		this.widthx = stack.getWidth();
		this.heighty = stack.getHeight();
		this.depthz = stack.getSize();

		// init other local variables according to the selected image
		labels3D = new int[widthx][heighty][depthz];
		distances3D = new float[widthx][heighty][depthz];
		

		this.gridSize = (grid < 5) ? 5 : grid;
		this.regul = (reg < 0) ? 0 : reg;

		// regul in range {0,1}
		this.factor = (regul * regul) * (float)(gridSize);	//TODO: check calculation


		estimTime = System.currentTimeMillis() - startTime;
		//printInfo("initializing variables took: " + estimTime/1000 + " s");
	}
	
	
	protected void convertImage() {
		
		long startTime, estimTime;
		startTime = System.currentTimeMillis();		
		
		switch (imagePlus.getType()) {		//FIXME

			case ImagePlus.COLOR_RGB:
				printInfo("converting RGB to cieLAB...");
				this.img3D = ConvertImage.rgb2cieLABfast(this.imagePlus);		
				
				estimTime = System.currentTimeMillis() - startTime;
				printInfo("converting RGB to cieLAB took: " + estimTime/1000 + " s");	
				
			break;
		
			// convert gray scale images

			case ImagePlus.GRAY8:
			case ImagePlus.GRAY16:
			case ImagePlus.GRAY32:

				printInfo("converting gray to cieLAB...");

				this.img3D = ConvertImage.gray2cieLAB3D(this.imagePlus);

				estimTime = System.currentTimeMillis() - startTime;
				printInfo("converting gray to cieLAB took: " + estimTime/1000 + " s");	
				
				break;
			default:
				printInfo("ERROR: Unsupported color space!");
				break;
		}
	}

	
	
	/**
	 * Initialize local variables and clusters. 
	 */	
	void initClusters() {
		printInfo("-> initializing clusters");
		
		// Calculate number of clusters
		// Math.ceil impacts a bigger gridSize more than a smaller gridSize
		int nbClusters = (int) (Math.ceil((float)widthx	 / (float)gridSize) * 
								Math.ceil((float)heighty / (float)gridSize) * 
								Math.ceil((float)depthz  / (float)gridSize));
		
		
		//System.out.println("Number of clusters: " + nbClusters);
		//clusterColor 
		clusterColor	= new int[nbClusters][img3D[0][0][0].length]; 		//img3D[0][0][0].length = 3 (nb of channels) 
		clusterPosition	= new int[nbClusters][3]; 	// 
		
		// Actual initialization of clusters
		
		int maxColumn = (int) Math.ceil(widthx  / (float)gridSize);	//maximum number of columns
		int maxRow	  = (int) Math.ceil(heighty / (float)gridSize);	//maximum number of rows 
		
		labels3D = new int [widthx][heighty][depthz];
		
		for (int x = 0; x < widthx; x++ ) {
			for (int y = 0; y < heighty; y++ ) {
				for (int z = 0; z < depthz; z++) {

					labels3D[x][y][z]= (int) ((x/gridSize) + (y/gridSize) * maxColumn + (z/gridSize) * maxColumn * maxRow);
				//	stack.setVoxel(x, y, z, labels3D[x][y][z]);
				}
			}
		}

		updateParallel();
		//clusterPosition = lowestGradient(image, clusterPosition);
		lowestGradient(clusterPosition);	
		//updateParallel(); //TODO: do i need this here again? 

		distGrid = null;



	}
	
	
	protected void lowestGradient (int[][] clusterPosition){
		printInfo("-> computing lowest gradient position...");
		
		long startTime, estimTime;
		startTime = System.currentTimeMillis();	
		
		
		Edges edges = new Edges();									// instantiate FeatureJ -> edges
		
		final Image image = Image.wrap(imagePlus);					// Change ImagePlus to Image to apply Edges()
		final double scale = 1;										// Scale of gausskernel (1st derivative)
		final boolean nonmaxsup = true;								// Boolean to determine if non local maxima gv should be set to 0
				
		Image gradientImage = edges.run(image, scale, nonmaxsup);	// Generate gradient Image
		ImagePlus gradientImagePlus = gradientImage.imageplus();	// Transform gradient Image back to an ImagePlus object
		
		double gv; 													// gv of clustercenter-voxel or connect 26 voxels around it
		double gvOld;												// variable to save gv of previous looping 
		int[][] connect26 = Connectivity3D.CONNECT26;			
		
		// Genereate new clusterPosition array with invariant clustercenters as reference
		int[][] clusterPositionNew = new int[clusterPosition.length][clusterPosition[0].length];
	
		for (int i = 0; i < clusterPosition.length; i++) {
			for (int j = 0; j < clusterPosition[0].length; j++) {
				clusterPositionNew[i][j] = clusterPosition[i][j];
			}
		}
			
		
		// Loop over all cluster-centers
		for (int c = 0; c < clusterPosition.length; c++) {			

			// make sure center is within distance 1 of image margin so that all 26 neighbor-voxels are within the image
			if((clusterPositionNew[c][0] * clusterPositionNew[c][1] * clusterPositionNew[c][2] == 0) 
				|| clusterPositionNew[c][0] >= (widthx-1) || clusterPositionNew[c][1] >= (heighty-1) || clusterPositionNew[c][2] >= (depthz-1)){ 
				
				continue;
		
			}


				// Save gv of cluster-center
				gv = gradientImagePlus.getStack().getVoxel(clusterPosition[c][0], clusterPosition[c][1], clusterPosition[c][2]);
				
				for (int n = 0; n < connect26.length; n++) {		// Loop over all neighbors of cluster-centers
					
					gvOld = gv;										// Save gv of cluster-center or previous loop-pass 
		 

					gv = gradientImagePlus.getStack().getVoxel(	clusterPositionNew[c][0] + connect26[n][0],  
																clusterPositionNew[c][1] + connect26[n][1], 
																clusterPositionNew[c][2] + connect26[n][2]);

					
					// if the new gv is darker ( = lower gradient) than the one from the previous loop pass, change clusterPosition coordinate to its voxel 
					if(gv < gvOld) {	
						clusterPosition[c][0] = clusterPositionNew[c][0] + connect26[n][0]; 
						clusterPosition[c][1] = clusterPositionNew[c][1] + connect26[n][1]; 
						clusterPosition[c][2] = clusterPositionNew[c][2] + connect26[n][2];
					}

				}	
				
				
		}	
		gradientImagePlus.close();
		
		this.clusterPosition = clusterPosition;
		
		estimTime = System.currentTimeMillis() - startTime;
		printInfo("computing lowest gradient position took: " + estimTime/1000 + " s");

	}
	

	
	/**
	 * Assign cluster index to each pixel in image according the given metric
	 */
	protected void updateParallel () {
		printInfo(" -> fast parallel update running...");

		long startTime, estimTime;
		startTime = System.currentTimeMillis();	
		
		// reset count of pixels that belong to a specific cluster. 
		int nbPixels[] = new int[clusterPosition.length];	// clusterPosition.length = Number of clusters. nbPixel[i] = nb of pixels in cluster i
		
		// reset all previous cluster centers
		for(int[] subarray : clusterColor) 		{	Arrays.fill(subarray, 0);	}
		for(int[] subarray : clusterPosition) 	{	Arrays.fill(subarray, 0);	}
		
		final ThreadUpdate[] threads = new ThreadUpdate[Threading.nbAvailableThread()];
		int deltaImg = (int) Math.ceil(widthx / (float)threads.length);		// divide image in slices (along x) with width deltaImg. -> assign as threads
		int endRange;

		for (int iThread = 0; iThread < threads.length; iThread++) {
			
			// Simultaneously run in as many threads as CPUs  
			threads[iThread] = new ThreadUpdate(img3D, gridSize, clusterPosition, clusterColor, labels3D, nbPixels); 
			// for all regular regions
			// because of a rounding the last has to cover rest of image

			endRange = (iThread + 1) * deltaImg;	// x-position at which the thread ends 
			threads[iThread].setRangeImg(iThread * deltaImg, endRange, 0, heighty, 0, depthz); // bX, eX, bY, eY, bZ, eZ; b = begin, e = end. 
		}
		

		Threading.startAndJoin(threads); 	//startAndJoin starts the passed thread, thus calls its run() method. here ThreadUpdate.run()
		
		int nb, cL, cA, cB, X, Y, Z;	
		
		// cycle over all clusters and divide them by nb pixels per cluster (= get mean position and color)
		for (int k = 0; k < clusterPosition.length; k++) {
			
			nb = 0; cL = 0; cA = 0; cB = 0; X = 0; Y = 0; Z = 0;
			
			// sum over threads
			for (int i = 0; i < threads.length; i++) {
				
				// returns count of pixels per cluster and thread as calculated by ThreadUpdate.run()
				// and sums it up over all threads, resulting in the total number of pixels per cluster....
				nb += threads[i].getNbPixels()[k];			
				
				
				//... same is done for each color and position 											
				// over all image channels
				cL += threads[i].getclusterColor()[k][0];	
				cA += threads[i].getclusterColor()[k][1];
				cB += threads[i].getclusterColor()[k][2];
				
				// over all positions
				X += threads[i].getClusterPositions()[k][0];
				Y += threads[i].getClusterPositions()[k][1];
				Z += threads[i].getClusterPositions()[k][2];

			}
			if (nb == 0) {		continue;	}
			nbPixels[k] = nb;
			
			// over all color channels
			clusterColor[k][0] = cL / nb;
			clusterColor[k][1] = cA / nb;
			clusterColor[k][2] = cB / nb;
			
			// over all positions
			clusterPosition[k][0] = X / nb;
			clusterPosition[k][1] = Y / nb;
			clusterPosition[k][2] = Z / nb;
		}		
		
		estimTime = System.currentTimeMillis() - startTime;
		//printInfo("parallel updating took: " + estimTime/1000 + " s");
	}
	
	
	/**
	 * compute distances of 2*gridSize search-Window from the cluster-center
	 */
	protected void computeDistGrid() {
		printInfo(" -> pre-computing the distance grid matrix...");

		long startTime, estimTime;
		startTime = System.currentTimeMillis();		

		int sz = 2 * gridSize + 1;				// search-window diameter
		distGrid = new float[sz * sz * sz]; 	// sz*sz*sz corresponds to the (2*gridSize)^3 search windows pixel-count
		float dx, dy, dz;
		
		// fill the array with distances to center
		for (int x = 0; x < sz; x++ ) {
			for (int y = 0; y < sz; y++ ) {
				for (int z = 0; z < sz; z++) {
					
					dx = x - gridSize + 1;
					dy = y - gridSize + 1;
					dz = z - gridSize + 1;

					distGrid[x * sz * sz + y * sz + z] = ((dx * dx) + (dy * dy) + (dz * dz)) * factor;	// iterate through every pixel in sz^3


				}
			}
		}

		estimTime = System.currentTimeMillis() - startTime;
		printInfo("computing dist grid took: " + estimTime/1000 + " s");
	}

	/**
	 * Count residual distance to nearest clusters by given metric
	 * 
	 * @return float returns a sum over all distances to nearest cluster
	 */
	protected float computeResidualError () {

		long startTime, estimTime;
		startTime = System.currentTimeMillis();	
		
		float err = 0;

		// cycle over all distances
		for (int x = 0; x < widthx; x++ ) {
			for (int y = 0; y < heighty; y++ ) {
				for (int z = 0; z < depthz; z++) {
					err += distances3D[x][y][z];
				}
			}
		}


		estimTime = System.currentTimeMillis() - startTime;
	//	printInfo("computing error took: " + estimTime/1000 + " s");
		return err;
	}
	
	
	/**
	 * Assign cluster index to each pixel in image according the given metric
	 */
	protected void assignmentParallel () {
		printInfo(" -> fast parallel assignement running...");
		long startTime, estimTime;
		startTime = System.currentTimeMillis();	
		
		// fill distances3D with highest possible positive Float value.  
		for(float[][] subarray : distances3D) {			
			for(float[] subsubarray: subarray) {
				Arrays.fill(subsubarray, Float.MAX_VALUE);	
			}
		}
		
		final ThreadAssignment[] threads = new ThreadAssignment[Threading.nbAvailableThread()];
		int deltaImg = (int) Math.ceil(widthx / (float)threads.length); 
		int endRange;
		
		for (int iThread = 0; iThread < threads.length; iThread++) {
			
			// Concurrently run in as many threads as CPUs  
			threads[iThread] = new ThreadAssignment(img3D, gridSize, distGrid, clusterPosition, clusterColor, distances3D, labels3D);
			// for all regular regions
			// because of a rounding the last has to cover rest of image

			endRange = (iThread + 1) * deltaImg;
			threads[iThread].setRangeImg(iThread * deltaImg, endRange, 0, heighty, 0, depthz); 
		
		}
		

		Threading.startAndJoin(threads); // startAndJoin calls ThreadAssignment.run()
		
		estimTime = System.currentTimeMillis() - startTime;
		printInfo("parallel assignment took: " + estimTime/1000 + " s");
	}
	

	
	
	/**
	 * Enforce Label Connectivity - Modified original code
	 * At the end of the clustering procedure, some orphaned pixels that do 
	 * not belong to the same connected component as their cluster center may 
	 * remain. To correct for this, such pixels are assigned the label of the 
	 * nearest cluster center using a connected components algorithm.
	 * 
	 * 1. finding an adjacent label for each new component at the start
	 * 2. if a certain component is too small, assigning the previously found
	 *    adjacent label to this component, and not incrementing the label.
	 */
	protected void enforceLabelConnectivity() {
		
		int maxlabel = 0;
		for (int i = 0; i < labels3D.length; i++) {
			for (int j = 0; j < labels3D[0].length; j++) {
				for (int k = 0; k < labels3D[0][0].length; k++) {
					
					if(maxlabel < labels3D[i][j][k]) {
						maxlabel = labels3D[i][j][k];
					}
				}
			}
		}

		System.out.println("maxLabel: "+ (maxlabel + 1));
	
		printInfo("maxLabel: "+ (maxlabel + 1));
		
		printInfo("Enforce label connectivity.");
		
//		ImageStack stack = imagePlus.getStack();
		long startTime, estimTime;
		startTime = System.currentTimeMillis();	
		
		// 6-connectivity	... TODO: check if this makes sense
		final int[] dx = {-1,  0,  0,  1,  0,  0 };	//x, y, x, y -> x, y, z, x, y, z
		final int[] dy = { 0, -1,  0,  0,  1,  0 };
		final int[] dz = { 0,  0, -1,  0,  0,  1 }; 	
		
//		final int dx[] = {-1,  0,  1,  0, -1,  1,  1, -1,  0, 0};
//		final int dy[] = { 0, -1,  0,  1, -1, -1,  1,  1,  0, 0};
//		final int dz[] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, 1};
		
		// 27-connectivity
//		final int dx[] = {0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0, 0, 0};
//		final int dy[] = {0, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 1, 1, 0,-1,-1,-1, 0, 1, 1, 0, 0};
//		final int dz[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0};

		
		
		// image size
		int imagevol = widthx * heighty * depthz;
		
		// area of initial superpixel
		int supvol = gridSize * gridSize * gridSize;	
		
		// create new array of labels and fill with -1
		int[][][] nlabels = new int[widthx][heighty][depthz];
		
		for(int[][] subarray : nlabels) { 		//TODO: doublecheck
			for(int[] subsubarray: subarray) {
				Arrays.fill(subsubarray, -1); 
			} 
		}

		// coordinates to run in the image
		int x, y, z;
		int lab = 0;
		int adjlabel = 0; //adjacent label
		// array of coordinates for all elements in the actual segment
		int[] xvec = new int[imagevol];
		int[] yvec = new int[imagevol];
		int[] zvec = new int[imagevol];
        int count;
		
        // cycle over all pixels in image
        for (int k = 0; k < depthz; k++) {
			for( int j = 0; j < heighty; j++ ) {
				for( int i = 0; i < widthx; i++ ) {
					
					if( nlabels[i][j][k] > -1) { 	continue; 	} // if this pixel already has been processed go to next looping
					
					nlabels[i][j][k] = lab;
					
					// Start a new segment
					xvec[0] = i;
					yvec[0] = j;
					zvec[0] = k;
					
					// find an adjacent label for later use if needed
					for( int n = 0; n < dx.length; n++ ) {
						
						x = xvec[0] + dx[n];
						y = yvec[0] + dy[n];
						z = zvec[0] + dz[n];
						
						if( (x >= 0 && x < widthx) && (y >= 0 && y < heighty) && (z >= 0 && z < depthz)) {	//checks if x,y,z is within image
							if(nlabels[x][y][z] >= 0) {
								adjlabel = nlabels[x][y][z];	//adjlabel: label of last 6-connected voxel >= 0
							}
					}
				}
			
				count = 1; // segment size
				// region growing method and storing pixels belongs to segment
				for( int c = 0; c < count; c++ ) {
					for( int n = 0; n < dx.length; n++ ) {
						x = xvec[c] + dx[n];
						y = yvec[c] + dy[n];
						z = zvec[c] + dz[n];
						// conditions if it is still the same segment
						if( (x >= 0 && x < widthx) && (y >= 0 && y < heighty) && (z >= 0 && z < depthz)){
							if( 0 > nlabels[x][y][z] && labels3D[i][j][k] == labels3D[x][y][z] ) {
								xvec[count] = x;
								yvec[count] = y;
								zvec[count] = z;
								nlabels[x][y][z] = lab;
								count++;
								
//								stack.setVoxel(x, y, z, 255);
//								try {
//								TimeUnit.MILLISECONDS.sleep(1);
//							} catch (InterruptedException e) {
//								// TODO Auto-generated catch block
//								e.printStackTrace();
//							}
//								imagePlus.setSlice(z);
//								imagePlus.updateAndDraw();
							}
						}
					}
				}
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				// shift by 2, which means that it reduces segments 4times smaller
				if(count <= supvol >> 2) {				//shift is calculated first, than relational operator.
					for( int c = 0; c < count; c++ ) {
						nlabels[xvec[c]][yvec[c]][zvec[c]] = adjlabel;
					}
					lab--;
				}
				lab++;
				}
			}
		}
		labels3D = nlabels;
		nbLabels = lab;
		//printInfo("number of labels: " + nbLabels);

		maxlabel = 0;
		for (int i = 0; i < labels3D.length; i++) {
			for (int j = 0; j < labels3D[0].length; j++) {
				for (int k = 0; k < labels3D[0][0].length; k++) {
					
					if(maxlabel < labels3D[i][j][k]) {
						maxlabel = labels3D[i][j][k];
					}
				}
			}
		}
//
//		
//		printInfo("maxLabel after enforcing labelconnectivity: "+ (maxlabel + 1));
		System.out.println("maxLabel after enforcing labelconnectivity: "+ (maxlabel + 1));
		
		
		
		estimTime = System.currentTimeMillis() - startTime;
		printInfo("enforcing label connectivity took: " + estimTime/1000 + " s");
	}
	
	
	/**
	 * used only for presenting the segmentation results
	 */
	protected void showSegmentation() {
		printInfo("jSLIC3D visualisation...");
		
		long startTime, estimTime;

		// show the ROI in ROI manager
		if (showROIs) {
			startTime = System.currentTimeMillis();	
			ij.IJ.log(" -> show ROI manager");
			//sp.getSegmentation().showOverlapROIs(imagePlus);
			
			getSegmentation().superPixelRois(labels3D);
			
			estimTime = System.currentTimeMillis() - startTime;
			printInfo("showing labels as color took: " + estimTime/1000 + " s");
		}
		
		
		// show the general Overlay
		if (clrOverlay != null) {
			startTime = System.currentTimeMillis();
			
			ij.IJ.log(" -> show contour overlap");

			getSegmentation().showOverlapContours(imagePlus, clrOverlay);
			estimTime = System.currentTimeMillis() - startTime;
			printInfo("showing contours took: " + estimTime/1000 + " s");
		}
		
	}

	

	/**
	 * gives segmentation with segmented indexes
	 * 
	 * @return int[Width][Height] returns indexes of segmented superpixels
	 */
	public Labelling3D getSegmentation() {
		return new Labelling3D(labels3D);
	}	
	
	
	/**
	 * gives the number of all various labels in segmentation, where the 
	 * max labels are {0,..,(n-1)}
	 * 
	 * @return int number of labels
	 */
	public int getNbLabels() {
		return this.nbLabels;
	}
	
	
	/**
	 * get the converted image in LAB colour space in case of RGB otherwise 
	 * only gray intensity values
	 * 
	 * @return int[Width][Height][channels]
	 */
	public int[][][][] getImage() {
		return this.img3D;
	}
	
	/**
	 * content of About frame
	 */
	public void showAbout() {
	    IJ.showMessage("About jSLIC...");	//TODO: add about information
	  } 
	
	/**
	 * 
	 * @param s string text
	 */
	private void printInfo(String s) {
		ij.IJ.log(s);
		ij.IJ.showStatus(s);
	}
	
	public static void main(final String... args) throws Exception{


		// Launch ImageJ
		final ImageJ ij = net.imagej.Main.launch(args);

		// open test image.
		IJ.run("T1 Head (2.4M, 16-bits)");		
		//IJ.run("/home/s321411/Bilder/test_stack.tif");

		
		// Launch "jSLIC_superpixels_3D"
		ij.command().run(Slic_ij_superpixels_3D.class, true);


	
	}
}
