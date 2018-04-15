/**
 * 
 */
package jSLIC3D;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
//import ij.gui.PolygonRoi;
//import ij.gui.Roi;
import ij.plugin.LutLoader;
//import ij.plugin.frame.RoiManager;
import ij.process.ColorProcessor;
//import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
//import ij.process.ShortProcessor;
import ij.process.StackConverter;
//import ij.process.StackProcessor;

import java.awt.Color;
//import java.io.FileNotFoundException;
//import java.io.PrintWriter;
//import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import jSLIC3D.Connectivity3D;

//import org.python.bouncycastle.crypto.generators.OpenSSLPBEParametersGenerator;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;


//import jSLIC3D.ConvertStructure;
//import mcib3d.image3d.*;
import mcib_plugins.tools.RoiManager3D_2;



/**
 * @class Labelling 2D
 * @version 0.1
 * @date 18/06/2013
 * @author Jirka Borovec <jiri.borovec@fel.cvut.cz>
 * @category image segmentation
 * 
 * @brief Derivation of an abstract class for Segmentation representation. 
 * This particular child handl only 2D images segmentations.
 * 
 */
public class Labelling3D /*extends Labelling*/ {
	int[] dims = null;
	// max label
	int maxLabel = 0;
	// histogram
	int[] hist = null;
	
	
	// labelling
	private int[][][] data = null;
	
	@Parameter
	LogService logService;
	
//	@Parameter
//	RoiManager3D_2 roiManager3D; 
	
	

//	/**
//	 * Construct empty labeling of given size w x h x d
//	 * 
//	 * @param w int width of new segmentation
//	 * @param h int height of new segmentation
//	 */
//	public Labelling3D(int w, int h, int d) {	
//		System.out.println("Labelling3D Constructor");
//		// rewrite data dimensions
//		dims = new int[3];
//		dims[0] = w;
//		dims[1] = h;
//		dims[2] = d;
//		
//		// init data array
//		data = new int[dims[0]][dims[1]][dims[2]];
//		for(int[][] subarray : data) {   
//			for(int[] subsubarray: subarray) {
//				Arrays.fill(subsubarray, 0);   
//			}
//		}
//	}
	
	/**
	 * Constructor
	 * 
	 * @param segm is new labeling matrix of int[width][height][depth]
	 */
	public Labelling3D(int[][][] segm) { //don't comment
		System.out.println("Labelling3D constructor 2");
		resetSegm(segm);
	}
	
	/**
	 * Reset the segmentation so that it copies new labeling and recomputes
	 * the histogram
	 * 
	 * @param segm is new labelling matrix of int[width][height][depth]
	 */
	public void resetSegm(int[][][] segm) {	//don't comment out
		System.out.println("Labelling3D resetSegm");
		copyData(segm);			
		computeHistogram();
	}

	/**
	 * Copy the data from input matrix to local representation and also update 
	 * the maximal label according new labeling
	 * 
	 * @param d is a new matrix of int[width][height][depth]
	 */
	protected void copyData(int[][][] d) {	//don't comment out
		System.out.println("Labelling3D copyData");
		// rewrite data dimensions
		dims = new int[3];
		dims[0] = d.length;
		dims[1] = d[0].length;
		dims[2] = d[0][0].length;	
		
		// init data array
		data = new int[dims[0]][dims[1]][dims[2]];
		
		// copy data
		for (int i = 0; i < d.length; i++) {
			for (int j = 0; j < d[i].length; j++) {
				for( int k = 0; k < d[i][j].length; k++){
					
					data[i][j][k] = d[i][j][k];
					
					if (d[i][j][k] > maxLabel) {
						maxLabel = d[i][j][k];
					}
				}
			}
		}
	}
	
//	/**
//	 * Gets the label in chosen position, in case out of segmentation it throws 
//	 * an exception
//	 * 
//	 * @param x int position in the first dimension
//	 * @param y int position in the second dimension
//	 * @param z int position in the third dimension
//	 * @return int label
//	 */
//	public int getLabel(int x, int y, int z) {
//		System.out.println("Labelling3D getLabel");
//		// check out of image	
//		if (x < 0 || y < 0 || z < 0 || x >= dims[0] || y >= dims[1] || z >= dims[2]) {
//			throw new IndexOutOfBoundsException();
//		}
//		return data[x][y][z];
//	}
	
//	/**
//	 * Sets the label in chosen position, in case out of segmentation it throws 
//	 * an exception
//	 * 
//	 * @param x int position in the first dimension
//	 * @param y int position in the second dimension
//	 * @param l int new label
//	 */
//	public void setLabel(int x, int y, int z, int l) {
//		System.out.println("Labelling3D setLabel");
//		maxLabel = -1;
//		// check out of image
//		if (x < 0 || y < 0 || z < 0 || x >= dims[0] || y >= dims[1] || z >= dims[2]) {
//			throw new IndexOutOfBoundsException();
//		}
//		// update histogram
//		if (hist != null) {
//			hist[ data[x][y][z] ] --;
//			hist[ l ] ++;
//		}
//		// assign
//		data[x][y][z] = l;
//		// update max label
//		if (l > maxLabel) {
//			maxLabel = l;
//		}
//	}
//	
//	/**
//	 * @see 
//	 */
//	//@Override
//	public int getLabel(int[] pos) {	
//		System.out.println("Labelling3D getLabel");
//		return getLabel(pos[0], pos[1], pos[2]);
//	}
//
//	/**
//	 * @see
//	 */
//	//@Override
//	public void setLabel(int[] pos, int l) {	
//		System.out.println("Labelling3D setLabel");
//		setLabel(pos[0], pos[1], pos[2], l);
//		
//	}
//	
//	/**
//	 * BE CAREFUL ABOUT THIS METHOD !!!
//	 * 
//	 * @return reference to data of the labelling
//	 */
//	public int[][][] getData() {
//		System.out.println("Labelling3D getData");
//		return data;
//	}
//
//	/**
//	 * 
//	 * @return
//	 */
//	public int getWidth() {
//		System.out.println("Labelling3D getWidth");
//		return data.length;
//	}
//
//	/**
//	 * 
//	 * @return
//	 */
//	public int getHeight() {
//		System.out.println("Labelling3D getHeight");
//		return data[0].length;
//	}
//	
//	/**
//	 * 
//	 * @return
//	 */
//	public int getDeptht() {	
//		System.out.println("Labelling3D getDepth");
//		return data[0][0].length;
//	}
//	
	

	/**
	 * @see sc.fiji.CMP_BIA.segmentation.structures.Labelling#computeHistogram()
	 */
	//@Override
	public int[] computeHistogram() {	//don't comment out
		System.out.println("Labelling3D computeHistogram");
		maxLabel = -1;
		// find new max labels 
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				for (int k = 0; k < data[i][j].length; k++) {
					if (maxLabel < data[i][j][k]) {
						maxLabel = data[i][j][k];
					}
				}
			}
		}

		// init hist. array		
		hist = new int[maxLabel+1];
		Arrays.fill(hist, 0);
		
		// compute histogram	
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				for (int k = 0; k < data[i][j].length; k++) {
					hist[ data[i][j][k] ] ++;
				}
			}
		}
		return hist;
	}
	
//	/**
//	 * @see sc.fiji.CMP_BIA.segmentation.structures.Labelling#showLabelling()
//	 */
//	//@Override
//	public void showLabelling() {
//		System.out.println("Labelling3D showLabelling");
//
//		
//		ImageStack segm2 = new ImageStack(dims[0], dims[1], dims[2]);	//TODO: make sure image does not exceed 16 bit (= short)
//		for (int x = 0; x < dims[0]; x++) {
//			for (int y = 0; y < dims[1]; y++) {
//				for (int z = 0; z < dims[2]; z++) {
//					segm2.setVoxel(x, y, z, data[x][y][z]);
//				}
//			}
//		}
//		
//		// create Processor to ImagePlus
//		ImagePlus img = new ImagePlus("Segmentation", segm2);
//		img.show();
//	}
//
//	/**
//	 * @see sc.fiji.CMP_BIA.segmentation.structures.Labelling#reLabel(int[] LUT)
//	 */
////	@Override
//	public void reLabel(int[] LUT) {
//		System.out.println("Labelling3D reLabel");
//		if ((maxLabel+1) != LUT.length) {
//			throw new IndexOutOfBoundsException("segmentation and new labelling LUT are not same.");
//		}
//		
//		assert LUT.length > 0;
//		maxLabel = LUT[0];
//		// find new max labels 
//		for (int i=1; i<LUT.length; i++) {
//			if (maxLabel < LUT[i]) {
//				maxLabel = LUT[i];
//			}
//		}
//		
//		// init hist. array		
//		hist = new int[maxLabel+1];
//		Arrays.fill(hist, 0);
//		
//		// relabel actual labelling and compute histogram	//TODO: done 1x
//		for (int i=0; i<data.length; i++) {
//			for (int j=0; j<data[i].length; j++) {
//				for (int k = 0; k < data[i][j].length; k++) {
//					data[i][j][k] = LUT[ data[i][j][k] ];
//					hist[ data[i][j][k] ] ++;
//				}
//			}
//		}
//	}
//	
//	/**
//	 * @see sc.fiji.CMP_BIA.segmentation.structures.Labelling#findSegmentsConnectivity(int[][])
//	 */
////	@Override
//	public int[][] findSegmentsConnectivity(int[][] neighbors) {
//		System.out.println("Labelling3D findSegmentsConnectivity");
//		return ConvertStructure.arrayLists2intMatrix( Connectivity3D.findSegmetNeighbors(data, maxLabel+1, neighbors) );
//	}
//
//	/**
//	 * @see sc.fiji.CMP_BIA.segmentation.structures.Labelling#findElementsBoundaries(int[][] neighborhood)
//	 */
//	//@Override
//	public ArrayList<ArrayList<int[]>> findElementsBoundaries(int[][] neighborhood) {
//		System.out.println("Labelling3D findElementsBoundaries");
//		return Connectivity3D.segmentBoundariesRaw(data, maxLabel+1, neighborhood);
//	}
//	
//	public ArrayList<ArrayList<int[]>> findElementsBoundariesPolygon() {
//		System.out.println("Labelling3D findElementsBoundariesPolygon()");
//		// treat all point and then simplify
//		//logService.info("   -> segment boundaries..");
//		ArrayList<ArrayList<int[]>> bounds = Connectivity3D.segmentBoundaries(data, maxLabel+1);
//		//logService.info("   -> simplify polygon...");
//		Connectivity3D.simplifyPolygon(bounds);
//		return bounds;
//	}
//
	/**
	 * @see sc.fiji.CMP_BIA.segmentation.structures.Labelling#showOverlapLabeling(ImagePlus img, float optically)
	 */
//	@Override
//	public void showOverlapLabeling(ImagePlus img, double optically) {
//		
//		if (img.getType() != ImagePlus.COLOR_RGB) {	//Add a warning if picture is not RGB
//
//			StackConverter stackConv = new StackConverter(img);
//			stackConv.convertToRGB();
//			
//			}
//		
//		System.out.println("Labelling3D showOverlapLabeling");
//		if ( ! checkImgAndSegmDims(img) ) {		return;		}
//				
//		// create LUT
//		Color clr = null;
//		Random rnd = new Random();
//		int[] lut = new int[maxLabel + 1];
//		int[][] lutRGB = new int[maxLabel + 1][3];
//		for (int i = 0; i <= maxLabel; i++) {
//			// segment colour in single integer
//			lut[i] = rnd.nextInt(255 * 255 * 255);
//			// decomposition by RGB components
//			clr = new Color(lut[i]);
//			lutRGB[i][0] = (int) (clr.getRed() * (1 - optically));
//			lutRGB[i][1] = (int) (clr.getGreen() * (1 - optically));
//			lutRGB[i][2] = (int) (clr.getBlue() * (1 - optically));
//		}	
//
//		// check if it is colour image
//		if (img.getType() != ImagePlus.COLOR_RGB) {		
//		//	logService.info("WARING: the image is not RGB image."); 
//			img.setProcessor( img.getProcessor().convertToRGB() );
//		} 
//
//		// pixel values (local)
//		int c[] = null; 	
//		// create colour segmentation
//		ImageStack stack = img.getStack();
//		ImageProcessor ip = img.getProcessor().convertToRGB();			//TODO: wrong processor?
//		ImageProcessor segm = new ColorProcessor(dims[0], dims[1]);
//
//		int maxSlices = img.getNSlices();
//		for (int z = 1; z <= maxSlices; z++) {
//			
//			ip = stack.getProcessor(z);
//			segm = stack.getProcessor(z);
//			
//			for (int i = 0; i < data.length; i++) {
//				for (int j = 0; j < data[i].length; j++) {
//					// segm.set(i, j, lut[ data[i][j] ]);
//					c = ip.getPixel(i, j, c); 
//					clr = new Color((int) (optically * c[0]) + lutRGB[data[i][j][z]][0],
//									(int) (optically * c[1]) + lutRGB[data[i][j][z]][1], 
//									(int) (optically * c[2]) + lutRGB[data[i][j][z]][2]);
//					segm.set(i, j, clr.getRGB() );
//					
//				}
//			}
//			//System.out.println(z);
//			stack.addSlice(segm);
//		}
//		
//
//		//stack.addSlice(segm);
//		//img.setStack(stack);
//		img.updateImage();
//		img.show();
//		
//	}

	/**
	 * @see sc.fiji.CMP_BIA.segmentation.structures.Labelling#showOverlapContours(ImagePlus img)
	 * 
	 * @param neighborhood is one of Connectivity3D.CONNECT4 or CONNECT8
	 * @example showOverlayContours(image, Connectivity3D.CONNECT4, Color.RED);
	 */
	//@Override
	public void showOverlapContours(ImagePlus img, java.awt.Color clr) {//don't comment out
		System.out.println("Labelling3D showOverlapContours");
		if ( ! checkImgAndSegmDims(img) ) {		return;		}	
		
		// check if it is a color image
		if (img.getType() != ImagePlus.COLOR_RGB) {	//Add a warning if picture is not RGB

			StackConverter stackConv = new StackConverter(img);
			stackConv.convertToRGB();
			
			}
		
		ImagePlus img2 = img.duplicate();
		
		ImageStack stack = img2.getImageStack();
		
		ArrayList<ArrayList<int[]>> coords = Connectivity3D.segmentBoundariesRaw(data, maxLabel+1, Connectivity3D.CONNECT26);

		for (int z = 1; z <= img2.getStackSize(); z++) {
			ImageProcessor sliceProcessor = stack.getProcessor(z);	
		//	System.out.println(z);
				
			// draw the contours
			for (int i = 0; i < coords.size(); i++) {	
				
				for (int j = 0; j < coords.get(i).size(); j++) {
					 if(coords.get(i).get(j)[2] == z) {
					
						
						sliceProcessor.set(coords.get(i).get(j)[0], 
											coords.get(i).get(j)[1], clr.getRGB());// clr.getRGB()); TODO: change to RGB
					 }
				}		
			}		
		}
		img2.updateAndRepaintWindow();	//This draws the pixels saved by ip.set in
		//img.updateAndDraw();
						
	//	img2.updateImage();
		img2.show();
	}
//
//	/**
//	 * 
//	 */
//	public void showOverlapROIs(ImagePlus img) {
//		System.out.println("Labelling3D showOverlapROIs");
////		System.out.println("labelling3D showOverlapROIs checkImgAndSegmDims " + checkImgAndSegmDims(img));
//		if ( ! checkImgAndSegmDims(img) ) {		return;		}
//		System.out.println("Labelling3D checkImgAnsSegmDims: everything okay with image dimensions");
//		// estimate boundaries
//		ArrayList<ArrayList<int[]>> coords = findElementsBoundariesPolygon();
//		
//		//int currentSlice = img.getCurrentSlice();
//		RoiManager manager = RoiManager.getInstance();
//		if (manager == null) {
//		    manager = new RoiManager();
//		}
//		
//		// for all boundaries
//		for (int i = 0; i < coords.size(); i++) {
//			// skip empty boundaries
//			if (coords.get(i) != null) {
//				FloatPolygon poly = new FloatPolygon();
//				for (int j = 0; j < coords.get(i).size(); j++) {
//					poly.addPoint( coords.get(i).get(j)[0] , coords.get(i).get(j)[1] );
//				}
//				PolygonRoi roi = new PolygonRoi(poly, Roi.POLYGON);
//				roi.setName("superpixel "+Integer.toString(i));
//				manager.addRoi(roi);
//				//rm.add(img, new PolygonRoi(poly, Roi.POLYGON), 0);
//				//manager.add(img, roi, currentSlice);
//			}
//		}
//					
//		img.updateAndDraw();
//	}
//	
	/**
	 * check dimensionality between image and labeling
	 * 
	 * @param img image to be compared with the labeling
	 * @return bool if the dimensions are consistent
	 */
	private boolean checkImgAndSegmDims(ImagePlus img) {//don't comment out
		System.out.println("Labelling3D checkImgAndSegmDims");
		
		if (dims[0] != img.getWidth() || dims[1] != img.getHeight() || dims[2] != img.getStackSize()) {
			//logService.info("ERROR: Inconsistent image and labeling size!");
			return false;
		} 
		return true;
	}

//	/**
//	 * Compute the overlap histogram of two segmentations
//	 * 
//	 * TODO: not optimise for large number of labels
//	 * 
//	 * @param lb is the other segmentation of the same dimension
//	 * @param shift is the relative shift of the second segmentation in relation to this one
//	 * @return int[this.maxLabel][lb.maxLabel] is sparse matrix
//	 */
//	public int[][] overlaps(Labelling3D lb, int[] shift) {
//		System.out.println("Labelling3D overlaps");
//		// inti the output array
//		int[][] overlap = new int[this.maxLabel+1][lb.maxLabel+1];
//
//		// variables depending on segmentation sizes
//		final int lDim = 3;	//TODO
//		int[] minDim = new int[lDim];
//		int[] lShiftA = new int[lDim];
//		int[] lShiftB = new int[lDim];
//		int[] end = new int[lDim];
//		
//		// do for both dimensions
//		for (int i=0; i<lDim; i++) {
//			// find minimal sizes of both segmentations
//			minDim[i] = (this.dims[i] < lb.dims[i]) ? this.dims[i] : lb.dims[i];
//			// find shifting for the second image
//			lShiftA[i] = (shift[i] >= 0) ? shift[i] : 0;
//			lShiftB[i] = (shift[i] < 0) ? -shift[i] : 0;
//			// find the ending of common range
//			end[i] = (shift[i] < 0) ? this.dims[i]+shift[i] : lb.dims[i]-shift[i];
//			// for case of overflow in the segm.
//			//end[i] = (end[i] > this.dims[i]) ? this.dims[i] : end[i]; 
//		}
//		
//		// go throw overlap of both segmentations
//		for (int i=0; i<end[0]; i++) {
//			for (int j=0; j<end[1]; j++) {
//				for (int k = 0; k < end[2]; k++) {
//					overlap[ this.data[i + lShiftA[0]] [j + lShiftA[1]] [k+lShiftA[2]]] [ lb.data[i+lShiftB[0]][j+lShiftB[1]] [k+lShiftB[2]]] ++;
//				}
//			}
//		}
//		
//		return overlap;
//	}

//	/**
//	 * @ see {@link sc.fiji.CMP_BIA.segmentation.structures.Labelling#findMultiClassBoundaryPoints(int[][])}
//	 */
////	@Override
//	public int[][] findMultiClassBoundaryPoints(int[][] neighbors) {
//		System.out.println("Labelling3D findMultiClassBoundaryPoints");
//		return ConvertStructure.arrayList2intMatrix( Connectivity3D.findBoundaryPoints(data, neighbors) );
//	}
//
//	/**
//	 * @ see {@link sc.fiji.CMP_BIA.segmentation.structures.Labelling#clone()}
//	 */
////	@Override
//	public Object clone() {
//		System.out.println("Labelling3D clone");
//		return new Labelling3D( this.data );
//	}

//	/**
//	 * 
//	 */
////	@Override	// TODO: 3D 
//	public void exportToFile(String path) {
//		PrintWriter out = null;
//		// create the string
//		String strDims = new String("Dims:");
//		for (int i = 0; i < dims.length; i++) {
//			strDims += " " + Integer.toString(dims[i]);
//		}
//		// IO process
//		try {
//			out = new PrintWriter(path, "UTF-8");
//			// write data
//			out.println(strDims);
//			for (int i=0; i<data.length; i++) {
//				for (int j=0; j<data[i].length; j++) {
//					out.print( Integer.toString( data[i][j] ) + " ");
//				}
//				out.println();
//			}
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		} catch (UnsupportedEncodingException e) {
//			e.printStackTrace();
//		} finally {
//			if (out != null) {
//				out.close();
//			}
//		}
//	}
//
//	/**
//	 * 
//	 */
////	@Override
//	public void printData() {
//		for (int i=0; i<data.length; i++) {
//			for (int j=0; j<data[i].length; j++) {
//				System.out.print( Integer.toString( data[i][j] ) + ", ");
//			}
//			System.out.println();
//		}
//	}
//	
//	
//	
	public void superPixelRois(int [][][] label3D) {
		System.out.println("Labelling3D superPixelRois");

		ImagePlus imagePlus;
		ImageStack imageStack;	
		LutLoader lutLoader = new LutLoader();
		int maxlabel;
		int bitdepth;
		long colorfactor;
		int[][][] labelVis = new int[label3D.length][label3D[0].length][label3D[0][0].length];
		
		// Find out biggest label and set bit-depth accordingly
		maxlabel = 0;
		
		for (int x = 0; x < label3D.length; x++) { 
			for (int y = 0; y < label3D[0].length; y++) { 
				for (int z = 0; z < label3D[0][0].length; z++) { 
					
					if (label3D[x][y][z] > maxlabel) {
						maxlabel = label3D[x][y][z];
					}
				}
			}
		}
		
		if (maxlabel < 256) {
			bitdepth = 8;
			colorfactor = 254;
		}
		else if(maxlabel < 65536) {	// 16 Bit
			bitdepth = 16;
			colorfactor = 65534;
		}
		else {						// 32 Bit
			bitdepth = 32;
			colorfactor = 4294967294l;
		}


			
			//normalize color for clearer visualization 
			for (int x = 0; x < label3D.length; x++) { 
				for (int y = 0; y < label3D[0].length; y++) { 
					for (int z = 0; z < label3D[0][0].length; z++) { 
						labelVis[x][y][z] = (int)((label3D[x][y][z]/(float)maxlabel) * colorfactor + 1);
					}
				}
			}



		
		imagePlus = IJ.createImage("Segmentation",label3D.length, label3D[0].length, label3D[0][0].length, bitdepth);
		imageStack = imagePlus.getStack();	
		
		// Visualize label3D in Image
		for (int x = 0; x < labelVis.length; x++) {
			for (int y = 0; y < labelVis[0].length; y++) {
				for (int z = 0; z < labelVis[0][0].length; z++) {

					imageStack.setVoxel(x, y, z, labelVis[x][y][z]);	
				}
			}
		}
		
		imagePlus.show();
		lutLoader.run("3-3-2 RGB");
		
//		Object[] args = {0};
//				
//		RoiManager3D_2 roiManager3D = new RoiManager3D_2();
//		
//		roiManager3D.handleExtension("Manager3D_AddImage", args);
	
		

	}
		
}
