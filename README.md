MULTISPECTRAL IMAGE PAN-SHARENING 

There are two source files in this repository. One is a Python (.py) source file, and the other is an IDL source file (.pro, the Interactive Data Language). I will give some details about the Python file first. The Python file pansharpen.py is a stand-alone program that is to be run on the UNIX command line environent. The first argument should be the string filename of a 1-band panchromatic Geotiff image file (.tif). The second argument should be a string filename of a 4-band multispectral geotiff image file (RGB,NIR). Both should have lowercase file extensions (.tif). The Python program does a number of tasks. If first uses GDAL tools to resample the multispectral Geotiff image file to the same higher dimensions as the panchromatic image Geotiff file using bicubic interpolation. This file is written to disk. 

Then, both the panchromatic and bicubic-resampled multispectral Geotiff image files are broken up into smaller Geotiff image tiles (using gdal_retile.py). These panchromatic-multispectral (bicubic resampled multispectral that is) tile pairs are then used to create pan-sharpened Geotiff tile files, one using the Brovey image sharpening technique, and the other using the FIHS sharpening algorithm (fast Intensity-Hue-Saturation method). Finally, these pan-sharpened tiles are put together into two large mosaic Geotiff image files (using gdal_merge.py). The entire tiling process is performed to conserve Python memory, and to avoid Python MemoryError's being thrown, as loading huge Numpy arrays into temporary memory in Python tends to flood memory. So the final outputs should be two pan-sharpened Geotiff image files (.tif): one using the Brovey algorithm, the other the FIHS. Each of of these Geotiff output files should have 4 layers (bands): RGB, and NIR. See below for an output example using a Landsat 8 scene. 

![Alt text](https://lh3.googleusercontent.com/-p8HiA4RuEJ8/VYmh_ttK-uI/AAAAAAAAADY/x3RXoBroTO8/w426-h390/pansharpeningExamples.jpg)

This Landsat scene was over western Greece. There, you can see the city of Argostoli on the island of Kefalonia. The upper-left panel is the original panchromatic 1-band grayscale image at 15m resolution. The upper-right panel is the original multispectral (RGB combination). The lower-left panel is the pan-sharpened RGB image in which the Brovey technique was used, and the lower-right panel the FIHS technique was used. 

usage: 
$ python pansharpen.py panfname.tif multispectral.tif 

The IDL program is an example of using IDL (with ENVI programming functions and utilities) to pansharpen a 4-band multispectral Geotiff image file (RGB,NIR) using the Gram-Schmidt pan-sharpening technique. To run that code, both IDL 8.0+ and ENVI 5.0+ should both be installed onto your UNIX/LINUX operating system. Please see the comments within the IDL source file itself for instructions on how to use the script. 

@author: 
Gerasimos Michalitsianos
Science Systems and Applications, Inc. 
June 2015 
