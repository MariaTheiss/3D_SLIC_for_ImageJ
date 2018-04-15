package jSLIC3D;

import java.util.Arrays;

import ij.Prefs;

public class Threading {


	/** 
	 * Create a Thread[] array as large as the number of processors available. 
     * From Stephan Preibisch's Multithreading.java class. See: 
     * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD 
     */  
    public static int nbAvailableThread() {  
    	return Prefs.getThreads();
        //return Runtime.getRuntime().availableProcessors();   
    } 
    
    /** 
     * Start all given threads and wait on each of them until all are done. 
     * From Stephan Preibisch's Multithreading.java class. See: 
     * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD 
     */  
    public static void startAndJoin(Thread[] threads) {  
        for (int ithread = 0; ithread < threads.length; ++ithread) {  
            threads[ithread].setPriority(Thread.NORM_PRIORITY);  
            threads[ithread].start();  
        }  
  
        try  {     
            for (int ithread = 0; ithread < threads.length; ++ithread)  
                threads[ithread].join();  
        } catch (InterruptedException ie) {  
            throw new RuntimeException(ie);  
        }  
    } 
	
}

/**
 * 
 * @author JB, MT
 *
 */
abstract class ThreadParticularImg3D extends Thread {

	protected int[][][][] 	img 				= null;		// source image
    protected int[][] 		clusterPositions	= null;		// cluster centers
    protected int[][] 		clusterColor 		= null;		// clusterColor
    protected int[][][] 	labels 				= null;		// labeling
    protected int 			gridSize;						// gridSize
    protected int[] 		rangeWidth, rangeHeight, rangeDepth;		// set range 
    	
    { setPriority(Thread.NORM_PRIORITY); }  
                	
    /**
     * initialization / copy reference to all needed variables 
     * 
     * @param im - image
     * @param cPos - cluster positions (cluster center positions)
     * @param cClr - cluster colors
     * @param lab - given labeling
     * @param gSize- grid-size
     */
    public ThreadParticularImg3D(int[][][][] im, int[][] cPos, int[][] cClr, int[][][] lab, final int gSize) { 
		img 				= im;			//Warum nicht .this? TODO
		clusterPositions 	= cPos;
		clusterColor 		= cClr;
		labels 				= lab;
    	gridSize 			= gSize;
	}
    
    /**
     * setting the particular cube in image to be processed
     * 
     * @param bW - begin in width dim
     * @param eW - end in width dim
     * @param bH - begin in height dim
     * @param eH - end in height dim
     * @param bD - begin in depth dim
     * @param eD - end in depth dim
     */
    public void setRangeImg(final int bW, final int eW, final int bH, final int eH, final int bD, final int eD) {
    	
    	rangeWidth		=	 new int[2];		//rangeWidth[0] contains x-start, rangeWidth[1] contains x-stop
    	rangeWidth[0] 	=	 (bW >= 0) ? bW : 0;
    	rangeWidth[1] 	=	 (eW < img.length) ? eW : img.length;
    	
    	rangeHeight 	=	 new int[2];		// rangeHeight  and rangeDepth go from 0 to end of image in this context
    	rangeHeight[0] 	=	 (bH >= 0) ? bH : 0;
    	rangeHeight[1] 	=	 (eH < img[0].length) ? eH : img[0].length;

    	rangeDepth		= 	 new int[2];
    	rangeDepth[0]	= 	 (bD >= 0) ? bD : 0;
    	rangeDepth[1]	= 	 (eD < img[0][0].length) ? eD : img[0][0].length;	
	}
}

/**
 * The particular thread for assignment in given region
 * @author JB
 *
 */
class ThreadAssignment extends ThreadParticularImg3D {  
    protected float[] 		distGrid	=	null; 	// precomputed distances
    protected float[][][] 	distances 	=	null;   // estimated distances
    		
    /**
     * initialization / copy reference to all needed variables 
     * 
     * @param im 	-  image
     * @param gSize -  grid-size
     * @param dGrid -  recomputed distance grid
     * @param cPos 	-  cluster-center positions
     * @param cClr 	-  cluster colors
     * @param dist 	-  the internal distances (whole image)
     * @param lab 	-  label of each pixel 	(whole image)
     */
    public ThreadAssignment(final int[][][][] im, final int gSize, final float[] dGrid, final int[][] cPos, final int[][] cClr, final float[][][] dist, int[][][] lab) {
		super(im, cPos, cClr, lab, gSize);
		distGrid	=	dGrid;
		distances	=	dist;
	}
    
    /**
     * the main body of the thread
     */
    @Override
    public void run() {  
    	// init
    	int 	xB, xE, yB, yE, zB, zE, i; // ?: xBeginning, xEnd, y", y", i: n pixels in sz*sz*sz 
		float 	dist, dL, dA, dB;
		int 	sz = 2 * gridSize + 1;
		
		// temporary variables - differences
		float distLAB;
		            	
    	// cycle over clusters and compute distances to all pixels in surrounding
		for (int k = 0; k < clusterPositions.length; k++) {

		
			xB = Math.max(0,				(int)(clusterPositions[k][0] - gridSize));	
			xE = Math.min(rangeWidth[1],	(int)(clusterPositions[k][0] + gridSize));	 
		
			yB = Math.max(0,				(int)(clusterPositions[k][1] - gridSize));
			yE = Math.min(rangeHeight[1],	(int)(clusterPositions[k][1] + gridSize));
			
			zB = Math.max(0, 				(int)(clusterPositions[k][2] - gridSize));
			zE = Math.min(rangeDepth[1], 	(int)(clusterPositions[k][2] + gridSize));
			
			// cycle over all pixels in 2 * gridSize region
			for (int x = xB; x < xE; x++ ) {				
				for (int y = yB; y < yE; y++) {
				
					i	 = 	(clusterPositions[k][0] - x		+ gridSize) * sz * sz;
					i 	+=	(clusterPositions[k][1] - y		+ gridSize) * sz;	 
					i 	+=	 clusterPositions[k][2]	- zB	+ gridSize;
					
					for (int z = zB; z < zE; z++, i-- ) {

						dL 		=	 img[x][y][z][0] - clusterColor[k][0];
						dA 		=	 img[x][y][z][1] - clusterColor[k][1];
						dB 		=	 img[x][y][z][2] - clusterColor[k][2];
						distLAB =	 (dL * dL) + (dA * dA) + (dB * dB);
						
						
						// Compute dist of each pixel
						dist = distLAB + distGrid[i];
						
			
						// if actual distance is smaller then the previous assign new label 
						if (dist < distances[x][y][z]) {	
							labels[x][y][z] 	= k;
							distances[x][y][z] = dist;
						}
					}
				}
			}			
		}
    	
    }
}

/**
 * The particular thread for update in given region
 * @author JB
 */
class ThreadUpdate extends ThreadParticularImg3D {

	protected int[] nbPixels = null;	// number of px per cluster
    protected int beginK, endK;   	 	// set range
	
	
    /**
     * initialization / copy reference to all needed variables 
     * 
     * @param im - image
     * @param cPos - clusters positions
     * @param cClr - cluster colors
     * @param lab - given labeling
     */
    public ThreadUpdate(final int[][][][] im, final int gSize, int[][] cPos, int[][] cClr, final int[][][] lab, int[] nbPx) {
		super(im, cPos, cClr, lab, gSize);	
    	nbPixels = nbPx;
	}
    
    public void setRange(final int start, final int stop) {
    	beginK 	= (start >= 0 ) ? start : 0;
    	endK	= (stop < nbPixels.length) ? stop : nbPixels.length;
	}

    
    @Override
    public void run() {
    	// Fill nbPixels with 0
		nbPixels = new int[nbPixels.length];
		Arrays.fill(nbPixels, 0);
		
		// Fill clusterColor with 0
		clusterColor = new int[nbPixels.length][3];
		for(int[] subarray : clusterColor)	 	{		Arrays.fill(subarray, 0);	}
		
		// fill clusterPositions with 0
		clusterPositions = new int[nbPixels.length][3];
		for(int[] subarray : clusterPositions) 	{		Arrays.fill(subarray, 0);	}
		
		int k;    	// cluster k
		
   		// cycle over all pixels in region
		for (int x = rangeWidth[0]; x < rangeWidth[1]; x++ ) {
			for (int y = rangeHeight[0]; y < rangeHeight[1]; y++) { 
				for (int z = rangeDepth[0]; z < rangeDepth[1]; z++) {
					
					k = labels[x][y][z];		// cluster-label of selected pixel (= to which cluster it belongs)
					
					// add up color - values for each cluster k
					clusterColor[k][0] += img[x][y][z][0];	// L	lightness
					clusterColor[k][1] += img[x][y][z][1];	// A	green-red
					clusterColor[k][2] += img[x][y][z][2];	// B	blue-yellow 
					
					
					// add up all positions in x, y and z for each cluster k 
					clusterPositions[k][0] += x;
					clusterPositions[k][1] += y;
					clusterPositions[k][2] += z;
					
					nbPixels[k] ++;		//count number of pixels 
				}
			}
		}	    	
    }
    
    public int[][] getclusterColor() {
		return clusterColor;
	}
    
    public int[][] getClusterPositions() {
		return clusterPositions;
	}
    
    public int[] getNbPixels() {
		return nbPixels;
	}
	
}