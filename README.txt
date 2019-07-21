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
INSTALLATION:
  $ git clone https://github.com/gerasimosmichalitsianos/pansharpen
  $ pip install pansharpen/
    or to upgrade:
  $ pip install pansharpen/ --upgrade
DIRECT USAGE:
  $ python pansharpen.py --panchromatic <PAN{.TIF}> --multispectral <MULTI{.TIF}>
    Options: 
      --version,       -v : display version help
      --help,          -h : display this usage messsage
      --panchromatic,  -p : pass in name of 1-band Geotiff holding 1-band panchromatic Geotiff image (high resolution)
      --multispectral, -m : pass in name of 3 or 4 band multispectral Geotiff image file (low-resolution)
GENERAL USAGE COMMAND LINE: 
     $ python-pansharpen 
      --panchromatic <PAN{.TIF}> 
      --multispectral <MULTI{.TIF}>
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
  These outputs should be in the same directory as the input files passed-in 
  via command-line.
AUTHOR: 
  Gersaimos A. Michalitsianos
  gerasimosmichalitsianos@gmail.com
  Last Updated: 21 July 2019

