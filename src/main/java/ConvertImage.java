/**
 * @file
 */
package jSLIC3D;

import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;

import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ChannelSplitter;
import ij.process.ImageProcessor;
import ij.process.StackConverter;

/**
 * @class Image Convertor
 * @version 0.1
 * @date 11/06/2013
 * @author Jirka Borovec <jiri.borovec@fel.cvut.cz>
 * @category image conversion
 * 
 * @brief converting an image into different colour spaces
 */
abstract public class ConvertImage {
	


//	/**
//	 * Convert whole image from RGB to LAB colour space
//	 * 
//	 * @param image is a ImageProcessor
//	 * @return int[width][height][3]
//	 */
//	public static int[][][] rgb2cieLAB (final ImageProcessor image) {
//		// check if it is RGB image
//		if (image.getNChannels() != 3) {
//			logService.info("Image is NOT RGB image, becase it has only "+ Integer.toString(image.getNChannels()) +" channels.");
//			return null;
//		}
//		
//		// create pixel buffer
//		int[][][] img = new int[image.getWidth()][image.getHeight()][3];
//		int c[] = null; // pixel values (local)
//		int lab[] = new int[3];
//		
//		for (int x=0; x<image.getWidth(); x++ ) {
//			for (int y=0; y<image.getHeight(); y++ ) {
//				c = image.getPixel(x, y, c); 
//				// Returns the pixel value at (x,y) as a 4 element array.
//				// RGB values are returned in the first 3 elements.
//		        ConvertColour.rgb2lab(c[0], c[1], c[2], lab);
//		        img[x][y][0] = lab[0];
//		        img[x][y][1] = lab[1];
//		        img[x][y][2] = lab[2];
//			}
//		}
//		
//		return img;
//	}
//
	/**
	 * Convert whole image from RGB to LAB colour space
	 * this fast version save already computed colour so each is computed only once
	 * 
	 * @param image is a ImageProcessor
	 * @return int[width][height][3]
	 */
	public static int[][][][] rgb2cieLABfast (final ij.ImagePlus image) {
		// check if it is an RGB image
		ImageProcessor processor = image.getProcessor();

		if (processor.getNChannels() != 3) {// imageProcessor returns 3 channels, when RGB, ImagePlus 1
			System.out.println("Image is NOT RGB, because it has only "+ Integer.toString(processor.getNChannels()) +" channels.");
			return null;
		}
		ImageStack stack = image.getStack();
		
		// saving already computed values
		int[][][][] LUT = new int[256][256][256][];
		
		int count = 0;
		
		// create pixel buffer
		int[][][][] img = new int[image.getWidth()][image.getHeight()][image.getNSlices()][3];
		//int c[] = null; // pixel values (local)
		int lab[] = new int[3];
	
		for (int z = 1; z <= image.getNSlices(); z++) {
			int c[] = null; // pixel values (local)
			
			processor = stack.getProcessor(z);	
				
			for (int x = 0; x < image.getWidth(); x++ ) {
				for (int y = 0; y < image.getHeight(); y++ ) {
	
					c = processor.getPixel(x, y, c); 
	
					// Returns the pixel value at (x,y) as a 4 element array.
					// RGB values are returned in the first 3 elements.
					if (LUT[c[0]][c[1]][c[2]] == null) {
						ConvertColour.rgb2lab(c[0], c[1], c[2], lab);
						LUT[c[0]][c[1]][c[2]] = lab.clone();
						count ++;
					} else {
						lab = LUT[c[0]][c[1]][c[2]].clone();
					}
			        img[x][y][z-1][0] = lab[0];
			        img[x][y][z-1][1] = lab[1];
			        img[x][y][z-1][2] = lab[2];
					}
				}
		}
		float rate = 	(float)count 
						/ (float)(image.getWidth()
						* image.getHeight()
						* image.getNSlices()) 
						* 100;
		System.out.println("   -> computed " + Integer.toString(count) + " colours = " + Float.toString(rate) + "%");
		return img;
	}
	
	/**
	 * Convert whole gray image own colour space
	 * 
	 * @param image is a ImageProcessor
	 * @return int[width][height][3]
	 */
	public static int[][][] gray2cieLAB(final ImageProcessor image) { //TODO: ImageProcessor
		// check if it is RGB image
		if (image.getNChannels() != 1) {
			System.out.println("Image is NOT gray scale, because it has only "+ Integer.toString(image.getNChannels()) +" channels.");
			return null;
		}
		
		// saving already computed values
		int[][] LUT = new int[256][];
				
		// create pixel buffer
		int[][][] img = new int[image.getWidth()][image.getHeight()][3];
		int[] vals = new int[4];
		int c; 						// pixel values (local)
		int lab[] = new int[3];		// int cRed, cGreen, cBlue;
		
		image.convertToByte(false);
		
		// over all pixels TODO: is it faster when i get Width and Height before entering the for-loop?
		for (int x = 0; x < image.getWidth(); x++ ) {
			for (int y = 0; y < image.getHeight(); y++ ) {

				c = image.getPixel(x, y, vals)[0];	// c ist der 8 bit grauwert des pixels
				
				if (LUT[c] == null) {
					ConvertColour.rgb2lab(c, c, c, lab);
					LUT[c] = lab.clone();
				} else {
					lab = LUT[c].clone();
				}
		        img[x][y][0] = lab[0];
		        img[x][y][1] = lab[1];
		        img[x][y][2] = lab[2];
				
			}
		}
		
		return img;
	}
	
	/**
	 *Trial to change color in 3D
	 */
	public static int[][][][] gray2cieLAB3D(final ij.ImagePlus image) {
		// check if it is RGB image
		if (image.getNChannels() != 1) {	// TODO: is a stack necessary here?
			System.out.println("Image is NOT gray scale, because it has only "+ Integer.toString(image.getNChannels()) +" channels.");
			return null;
		}
		StackConverter stackConv = new StackConverter(image);
		stackConv.convertToGray8();	// convert image to 8 bit grey-scale

		ImageStack stack = image.getStack();		
		
		// saving already computed values
		int[][] LUT = new int[256][];


		// get Stack-dimensions
		int stackWidthx = stack.getWidth();
		int stackHeighty = stack.getHeight();
		int stackDepthz = stack.getSize();

		// create pixel buffer
		int[][][][] img = new int[stackWidthx][stackHeighty][stackDepthz][3];
		double gv_double;	// gv as double value
		int c; 				// pixel values (local)
		int lab[] = new int[3];		// int cRed, cGreen, cBlue;

		// over all pixels 
		for (int x = 0; x < stackWidthx; x++ ) {
			for (int y = 0; y < stackHeighty; y++ ) {
				for (int z = 0; z < stackDepthz; z++) {
	
					
					gv_double = stack.getVoxel(x, y, z);
					c = (int)gv_double;
					
					
					if (LUT[c] == null) {
						ConvertColour.rgb2lab(c, c, c, lab);
						LUT[c] = lab.clone();
					} else {
						lab = LUT[c].clone();
					}
			        img[x][y][z][0] = lab[0];
			        img[x][y][z][1] = lab[1];
			        img[x][y][z][2] = lab[2];
				}
			}
		}
		System.out.println("ConvertImage: end grey2cielab3d");
		return img;
	}
	
	
	
	/**
	 * Convert whole gray image own color space
	 * 
	 * @param image is a ImageProcessor
	 * @return int[width][height][1]
	 */
	public static int[][][] gray2bright(final ImageProcessor image) {	//TODO: ImageProcessor
		// check if it is RGB image
		if (image.getNChannels() != 1) {
			System.out.println("Image is NOT gray image, becase it has only "+ Integer.toString(image.getNChannels()) +" channels.");
			return null;
		}
				
		// create pixel buffer
		int[][][] img = new int[image.getWidth()][image.getHeight()][1];
		int[] vals = new int[4];
		
		image.convertToByte(false);

		// over all pixels
		for (int x=0; x<image.getWidth(); x++ ) {
			for (int y=0; y<image.getHeight(); y++ ) {
				img[x][y][0] = image.getPixel(x, y, vals)[0];
			}
		}
		
		return img;
	}
	
//	
//	/**
//	 * Convert whole RGB image own colour space by brightness
//	 * 
//	 * @param image is a ImageProcessor
//	 * @return float[width][height]
//	 */
//	public static float[][] rgb2bright(final ImageProcessor image) {
//		// check if it is RGB image
//		if (image.getNChannels() != 3) {
//			System.out.println("Image is NOT RGB image, becase it has only "+ Integer.toString(image.getNChannels()) +" channels.");
//			return null;
//		}
//		
//		// create pixel buffer
//		float[][] img = new float[image.getWidth()][image.getHeight()];		
//		int c[] = null; // pixel values (local)
//		
//		// cycle over whole image and by labels add current value to given cluster center
//		for (int x=0; x<image.getWidth(); x++ ) {
//			for (int y=0; y<image.getHeight(); y++ ) {
//				c = image.getPixel(x, y, c); 
//				img[x][y] = ConvertColour.rgb2bright(c[0], c[1], c[2]);
//			}
//		}
//		
//		return img;
//	}
//
	
	
	
}
