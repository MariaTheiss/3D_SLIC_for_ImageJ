# 3D SLIC for ImageJ
3D version of jSLIC: 
https://imagej.net/plugins/cmp-bia-tools/ 
For the original implementation, see also Borovec and Kybic, 2014: 
https://core.ac.uk/download/pdf/47181098.pdf

I finished this work in 2018 and it did generate 3D superpixels as a Fiji plugin. 
There were remaining issues with anisotropic pixel-sizes. 
This code might not function out of the box anymore. 

## Overview of .java files: 

### Slic_ij_superpixels_3D.java 
&rarr; main class

### Threading.java
&rarr; contains functionality for threading, including assignment and update.

### ConvertImage.java
&rarr; coordinates color conversion. Determines image color space and calls color conversion functions from ConvertColour.java.

### ConvertColour.java
&rarr; converts between color spaces on a pixel level. 

### Connectivity3D.java
&rarr; includes functionality for connectivity.

### Labelling3D.java
&rarr; visualization functions.
