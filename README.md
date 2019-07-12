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
       This file is written to disk. For this resampled multispectral Geotiff image
       file and the panchromatic Geotiff image file, the following pan-sharpening
       algorithms are applied:
         (1) Brovey, 
         (2) Principal Component Analysis (PCA),
         (3) FIHS (Fast Intensity Hue Saturation),
         (4) Wavelet
       
       
###### usage message:
       
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
          
          NOTE: These outputs should be in the same directory as the input files passed-in 
          via command-line.
       PYTHON VERSION:
         Supports Python 2.7.x
      AUTHOR: 
        Gerasimos A. Michalitsianos
        gerasimosmichalitsianos@gmail.com
        July 2019
        
###### Sample Outputs
        
        The below sample images show the results of this algorithm using Landsat 8 imagery over Skala (Σκάλα), 
        Greece. Skala is a small Greek town found on Greece's island of Kefalonia in western Greece.

![Alt text](https://i.imgur.com/QYxruGN.png)

       Left: original panchromatic image.
       Center: pan-sharpened RGB image using Brovey technique.
       Right: pan-sharpened RGB image using FIHS (Fast Intensity Hue Saturation) technique.

![Alt text](https://i.imgur.com/CUJt4JK.png)

       Left: original RGB low-resolution image.
       Center: pan-sharpened RGB image using Wavelet technique.
       Right: pan-sharpened RGB image using PCA (Principal Component Analysis) technique.

###### usage: 
       $ python pansharpen.py --panchromatic <PAN{.TIF}> --multispectral <MULTI{.TIF}>

###### @author: 
       Gerasimos Michalitsianos
       gerasimosmichalitsianos@gmail.com
       July 2019
