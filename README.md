###### MULTISPECTRAL IMAGE PAN-SHARENING 

       There is one source file in this repository: a single Python (.py) 
       source file. The Python file pansharpen.py is a stand-alone program that 
       is to be run on the UNIX/Linux command line environent. The first required 
       argument should be the filename of a 1-band panchromatic Geotiff image file
       (.tif extension). The second argument should be a string filename of a 3 or 
       4-band multispectral geotiff image file (RGB,NIR bands, in that order).
       The Python program does a number of tasks. If first uses GDAL tools to 
       resample the multispectral Geotiff image file to the same higher dimensions 
       as the panchromatic image Geotiff file using bicubic interpolation. 
       This file is written to disk. Then, both the resampled multispectral and 
       panchromatic Geotiffs are sliced-up into tiles (using gdal_retile.py). For 
       each pair, smaller pan-sharpened multispectral Geotiffs are created using 
       4 methods of pan-sharpening (resulting in 4 new MS Geotiffs for each tile): 
         (1) Brovey, 
         (2) Principal Component Analysis (PCA),
         (3) FIHS (Fast Intensity Hue Saturation),
         (4) Wavelet
       Later, for each method of pan-sharpening, all tiles are merged together (mosaic)
       to form 4 lare mosaics (using gdal_merge.py) for the 4 methods of pan-sharpening.
       
       NAME: 
         pansharpen.py
       DESCRIPTION:
        This program performs pansharpening of satellite imagery. It is meant to be 
        run on the command-line on UNIX-like operating systems. The two primary inputs
        are (1) a 3 or 4 band multispectral geotiff containing the red, green, blue, 
        and NIR bands (NIR is optional) and (2) a 1-band geotiff containing higher-
        resolution greyscale panchromatic image data. It is assumed that both of 
        these two Geotiff inputs are "clipped" to the same rectangular geographic 
        bounding-box. Four methods of pan-sharpening are used: Brovey, Fast Intensity
        Hue Saturation (FIHS), Wavelet, and Principal Component Analysis (PCA).
      USAGE:
        $ python pansharpen.py --panchromatic <PAN{.TIF}> --multispectral <MULTI{.TIF}>
         Options: 
          --version,       -v : display version help
          --help,          -h : display this usage messsage
          --panchromatic,  -p : pass in name of 1-band Geotiff holding 1-band panchromatic Geotiff image (high resolution)
          --multispectral, -m : pass in name of 3 or 4 band multispectral Geotiff image file (low-resolution)
          --gdalretile,    -r : path to gdal_retile.py (GDAL)
          --gdalmerge,     -s : path to gdal_merge.py (GDAL)
      EXAMPLE USAGE:
        This program can be run at the Linux/UNIX command-line.
        $ multispectralGeotiff=LC08_L1TP_185033_20170712_20170726_01_T1_MULTI_TOA_3BAND.TIF
        $ panchromaticGeotiff=LC08_L1TP_185033_20170712_20170726_01_T1_B8_TOA.TIF
        $ python pansharpen.py --panchromatic $panchromaticGeotiff --multispectral $multispectralGeotiff
      OUTPUTS: 
        When the program is complete, there should be 4 new Geotif image files: 
          (1) a 3 or 4 band multispectral Geotiff image created using Brovey pan-sharpening
          (2) a 3 or 4 band multispectral Geotiff image created using FIHS pan-sharpening
          (3) a 3 or 4 band multispectral Geotiff image created using Wavlet pan-sharpening
          (4) a 3 or 4 band multispectral Geotiff image created using PCA pan-sharpening
          These outputs should be in the same directory as the input files passed-in 
          via command-line.
      AUTHOR: 
        Gersaimos A. Michalitsianos
        gerasimosmichalitsianos@gmail.com
        February 2019

![Alt text](https://lh3.googleusercontent.com/-p8HiA4RuEJ8/VYmh_ttK-uI/AAAAAAAAADc/1220R530qfM/w800-h800/pansharpeningExamples.jpg)

       This Landsat scene was over western Greece. There, you can see the city of Argostoli on 
       the island of Kefalonia. The upper-left panel is the original panchromatic 1-band grayscale 
       image at 15m resolution. The upper-right panel is the original multispectral (RGB combination).
       The lower-left panel is the pan-sharpened RGB image in which the Brovey technique was used,
       and the lower-right panel the FIHS technique was used. 

###### usage: 
       $ python pansharpen.py --panchromatic <PAN{.TIF}> --multispectral <MULTI{.TIF}>

###### @author: 
       Gerasimos Michalitsianos
       SGT., Inc. 
       February 2018 
