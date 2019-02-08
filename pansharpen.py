#!/usr/bin/env python
import os
import sys
import gc 
import numpy as np
import warnings as warn 
import getopt
import subprocess
import pywt
from osgeo import osr, gdal, gdalconst

class pansharpen(object): 

  def __init__(self,pan,multi):
    '''function __init__(self,pan,multi):
    This is a constructor method to create an object of this class.
    Instantiate an object of this class with the full path filename
    strings for a 1-band panchromatic Geotiff image file and its 
    corresponding 3 or 4-band multispectral Geotiff image file 
    (with bands Red,Green,Blue, and NIR (optional) in that order.

    Args:
      pan (str): Name of panchromatic Geotiff image file.
      multi (str): Name of multispectral Geotiff image file.
    Returns: 
      pancharpen: an object of the pansharpen class.

    '''
    self.pan = pan
    self.multi = multi

  def pansharpenPCA(self): 
    '''function pansharpenPCA(self):
    This is an instance method that returns a Python list of 3 or 4
    NumPy arrays containing the pan-sharpened Red,Green,Blue, and 
    optionally, NIR bands. These bands will have been created using 
    the principal component analysis (PCA) pan-sharpening.
    
    Returns: 
      list: Python list[] containing 3 or 4 NumPy arrays using PCA method.
    '''

    # read Panchromatic,Multispectral Geotiffs into GDAL objects
    dsPan   = gdal.Open( self.pan )
    dsMulti = gdal.Open( self.multi )

    # get the number of rows (Y) and number of columns (X)
    nrows,ncols = dsPan.RasterYSize , dsPan.RasterXSize

    # read Panchromatic,Red,Green,Blue bands into 2D arrays
    pan     =   dsPan.GetRasterBand(1).ReadAsArray().astype(float)
    blue    = dsMulti.GetRasterBand(3).ReadAsArray().astype(float)
    green   = dsMulti.GetRasterBand(2).ReadAsArray().astype(float)
    red     = dsMulti.GetRasterBand(1).ReadAsArray().astype(float)

    nanlocs = np.where( pan == 0.0 )

    # create multiband image object (can be either 3 or 4 bands)
    if dsMulti.RasterCount == 4:

      NIR = dsMulti.GetRasterBand(4).ReadAsArray().astype(float)
      image = np.zeros( (nrows,ncols,4),dtype=np.float32)
      image[:,:,0] = red
      image[:,:,1] = green
      image[:,:,2] = blue 
      image[:,:,3] = NIR

    else: 
      image = np.zeros( (nrows,ncols,3) , dtype=np.float32)
      image[:,:,0] = red
      image[:,:,1] = green 
      image[:,:,2] = blue

    # perform PCA
    m,n,d=image.shape
    M = np.reshape( image, (m*n,d))
    [ PCAData , PCAMap ] = pansharpen.multiPCA( M,d )
    
    if PCAData.size == 1: 
      PCAData = np.zeros( (m,n,d)) + np.array(PCAData)[0,0]
      PCAMap['M'] = PCAMap['M'].T
    else: 
      PCAData = np.reshape( PCAData , (m,n,d))

    nanlocs=np.where(pan==0.0)
    pan = np.where( pan == 0.0 , np.nan, pan)
    pan = np.ma.array( pan, mask=np.isnan(pan))

    F = PCAData
    Fmasked=np.zeros(F.shape)
    for band in range(F.shape[2]):
      Fmasked[:,:,band] = F[:,:,band]
      Fmasked[:,:,band][nanlocs] = 0.0

    Fmasked = np.where(Fmasked == 0.0 , np.nan , Fmasked)
    Fmasked = np.ma.array( Fmasked, mask=np.isnan(Fmasked))

    P = ( pan - pan.mean() ) * ( Fmasked / pan.std()).std() - Fmasked.mean()
    F[:,:,0] = P
    F = pansharpen.inversePCA( PCAMap['M'] , np.reshape( F, (m*n,d)), PCAMap['mean'])
    sharpened = np.reshape(F,(m,n,d))

    outred   = sharpened[ :,:,0 ] 
    outgreen = sharpened[ :,:,1 ]
    outblue  = sharpened[ :,:,2 ] 

    outred[nanlocs]   = 0.0
    outgreen[nanlocs] = 0.0
    outblue[nanlocs]  = 0.0

    if dsMulti.RasterCount == 4: 
      outNIR = sharpened[ :,:,3 ]
      outNIR[nanlocs] = 0.0
      return [ outred, outgreen, outblue, outNIR ]
    elif dsMulti.RasterCount == 3:
      return [ outred, outgreen, outblue ] 
    else: 
      return []

  @staticmethod
  def inversePCA( E,P,MeanV ):
    m,n = P.shape
    dim = MeanV.shape[0]
    M   = np.zeros((m,dim))
    
    for k in range(0,dim):  M[ 0:m, k  ] = MeanV[k]
    return np.ma.dot( P, E.T ) + M

  @staticmethod
  def multiPCA( X , no_dims=None ):
    '''
    function multiPCA( X, no_dims=None): 
    This function does necessary computations as part of the 
    principal component analysis (PCA) pan-sharpening method.

    Args:
      X (numpy.ndarray): 3 or 4-layer data cube (3D numpy array).
      no_dims (int): number of dimensions. Default is 2. 
    '''

    # create output dictionary to hold PCA results
    X = np.where( X == 0.0 , np.nan, X)
    X = np.ma.array( X , mask=np.isnan(X))
    mapping = {}

    # establish dims. if not passed-in as argument
    if no_dims is None: 
      no_dims = 2
    
    mapping['mean'] = np.ma.mean( X , axis=0)
    X = X - np.tile( mapping['mean'] , (X.shape[0] , 1 )   )

    # compute covariance matrix
    if X.shape[1] < X.shape[0]:
      C = np.ma.cov( X , rowvar=False)
    else:
      C = ( 1 / X.shape[0] ) * np.dot( X , X.T )

    # perform eigenvalue decomposition of C
    C[np.isnan(C)] = 0
    C[np.isinf(C)] = 0

    eigenvalues,M = np.linalg.eigh( C, UPLO='U')
    ind     = np.arange( 0 , eigenvalues.shape[0] , 1)
    ind     = np.flip(ind)

    if no_dims > M.shape[1]: # second axis
      no_dims =  M.shape[1]

    M = M[:,ind[0:no_dims]]
    eigenvalues = eigenvalues[0:no_dims]

    if not ( X.shape[1] < X.shape[0] ):
      with warn.catch_warnings():
        warn.filterwarnings('ignore',category=RuntimeWarning)  
        M = np.dot(X.T,M) * np.tile(  (1.0 / np.sqrt( X.shape[0] * eigenvalues )).T , ( X.shape[1], 1))

    mappedX  = np.ma.dot( X , M )
    mapping['M']       = M
    mapping['lambda']  = eigenvalues
    return [ mappedX, mapping ]

  def pansharpenFIHS(self):
    '''function pansharpenFIHS(self):
    This is an instance method that returns a Python list of 3 or 4
    NumPy arrays containing the pan-sharpened Red,Green,Blue, and 
    optionally, NIR bands. These bands will have been created using 
    the Fast Intensity Hue Saturation (FIHS) pan-sharpening method
    
    Returns: 
      list: Python list[] containing 3 or 4 NumPy arrays using FIHS method.
    '''
    # read Panchromatic,Multispectral Geotiffs into GDAL objects
    dsPan   = gdal.Open(self.pan)
    dsMulti = gdal.Open(self.multi)

    # read the Panchromatic,Red,Green,Blue bands into 2D arrays
    pan     =   dsPan.GetRasterBand(1).ReadAsArray().astype(float)
    blue    = dsMulti.GetRasterBand(3).ReadAsArray().astype(float)
    green   = dsMulti.GetRasterBand(2).ReadAsArray().astype(float)
    red     = dsMulti.GetRasterBand(1).ReadAsArray().astype(float)

    if dsMulti.RasterCount == 3:

      L = ( red+green+blue ) / 3.0
      redsharp   = red   + ( pan   - L )
      greensharp = green + ( green - L )
      bluesharp  = blue  + ( blue  - L )
      dsPan,dsMulti=None,None
      return [redsharp,greensharp,bluesharp]

    elif dsMulti.RasterCount == 4: 

      NIR   = dsMulti.GetRasterBand(4).ReadAsArray().astype(float)
      L = ( red+green+blue+NIR ) / 4.0

      redsharp   = red   + ( pan   - L )
      greensharp = green + ( green - L )
      bluesharp  = blue  + ( blue  - L )
      NIRsharp   = NIR   + ( NIR   - L )

      dsPan,dsMulti=None,None
      return [redsharp,greensharp,bluesharp,NIRsharp]

    else: 
      dsPan,dsMulti = None,None
      return []

  def pansharpenWavelet(self):
    '''function pansharpenWavelet(self):
    This is an instance method that returns a Python list of 3 or 4
    NumPy arrays containing the pan-sharpened Red,Green,Blue, and 
    optionally, NIR bands. These bands will have been created using 
    the Wavelet  pan-sharpening method
    
    Returns: 
      list: Python list[] containing 3 or 4 NumPy arrays using wavelet method.
    '''
    # read Panchromatic,Multispectral Geotiffs into GDAL objects
    dsPan = gdal.Open(self.pan)
    dsMulti = gdal.Open(self.multi)

    # read Panchromatic,Red,Green,Blue bands into 2D NumPy arrays
    pan  = dsPan.GetRasterBand(1).ReadAsArray()
    blue    = dsMulti.GetRasterBand(3).ReadAsArray().astype(float)
    green   = dsMulti.GetRasterBand(2).ReadAsArray().astype(float)
    red     = dsMulti.GetRasterBand(1).ReadAsArray().astype(float)
    d = dsMulti.RasterCount
    nrows,ncols = pan.shape
    
    if d == 3: 
      
      image = np.zeros((nrows,ncols,d),dtype=np.float32)
      image[:,:,0] = red
      image[:,:,1] = green
      image[:,:,2] = blue

    elif d == 4: 
      NIR = dsMulti.GetRasterBand(1).ReadAsArray().astype(float)
      image[:,:,3] = NIR
    else: 
      dsPan,dsMulti = None,None
      return []

    level = 0
    wavelet_type = 'haar'

    coeffs = pywt.wavedec2( pan, wavelet=wavelet_type, level=level)
    panvec,coeff_slices,coeff_shapes = pywt.ravel_coeffs(coeffs)
    reconstvec = np.tile( panvec.T , (d,1)).T

    n=panvec.shape[0]
    lowresvec = np.zeros((n,d),dtype=np.float32)

    for band in range(d):
      lowresCoeffs = pywt.wavedec2( image[:,:,band], wavelet=wavelet_type, level=level)
      lowresArr,arrSlices = pywt.coeffs_to_array(lowresCoeffs)
      lowresvec[:,band] = np.reshape(lowresArr,(nrows*ncols,))

    for j in range( 0 , coeff_shapes[0][0] * coeff_shapes[0][1] ):
      reconstvec[ j,:] = lowresvec[j,:]

    sharpened = np.zeros((nrows,ncols,d),dtype=np.float32)
    for band in range(d):
      p = np.reshape( reconstvec[:,band], (nrows,ncols))
      fcoeffs = pywt.wavedec2(p,wavelet_type,level=level)
      out=pywt.waverec2(fcoeffs,wavelet_type)
      sharpened[:,:,band] = out

    redsharp = sharpened[:,:,0]
    greensharp = sharpened[:,:,1]
    bluesharp = sharpened[:,:,2]
    
    if d == 4:
      NIRsharp = sharpened[:,:,3]
      return [redsharp,greensharp,bluesharp,NIRsharp]
    elif d == 3: 
      return [redsharp,greensharp,bluesharp]
    else: 
      return []

  def pansharpenBrovey(self):
    '''function pansharpenWavelet(self):
    This is an instance method that returns a Python list of 3 or 4
    NumPy arrays containing the pan-sharpened Red,Green,Blue, and 
    optionally, NIR bands. These bands will have been created using 
    the Brovey pan-sharpening method
    
    Returns: 
      list: Python list[] containing 3 or 4 NumPy arrays using Brovey method.
    '''
    # read Panchromatic,Multispectral Geotiffs into GDAL objects
    dsPan   = gdal.Open(self.pan)
    dsMulti = gdal.Open(self.multi)

    # read the Panchromatic,Red,Green,Blue bands into 2D arrays
    pan     =   dsPan.GetRasterBand(1).ReadAsArray().astype(float)
    blue    = dsMulti.GetRasterBand(3).ReadAsArray().astype(float)
    green   = dsMulti.GetRasterBand(2).ReadAsArray().astype(float)
    red     = dsMulti.GetRasterBand(1).ReadAsArray().astype(float)

    with warn.catch_warnings():
      warn.filterwarnings('ignore',category=RuntimeWarning)

      if dsMulti.RasterCount == 3: 

        redsharp   = np.multiply( np.true_divide( red, red+green+blue)   , pan )
        greensharp = np.multiply( np.true_divide( green, red+green+blue) , pan )
        bluesharp  = np.multiply( np.true_divide( blue, red+green+blue)  , pan )
        dsPan,dsMulti=None,None
        return [redsharp,greensharp,bluesharp]

      elif dsMulti.RasterCount == 4:

        NIR   = dsMulti.GetRasterBand(4).ReadAsArray().astype(float)
        redsharp   = np.multiply( np.true_divide( red, red+green+blue+NIR)   , pan )
        greensharp = np.multiply( np.true_divide( green, red+green+blue+NIR) , pan )
        bluesharp  = np.multiply( np.true_divide( blue, red+green+blue+NIR)  , pan )
        NIRsharp   = np.multiply( np.true_divide( NIR, red+green+blue+NIR)   , pan )
      
        dsPan,dsMulti=None,None
        return [redsharp,greensharp,bluesharp,NIRsharp]

      else:
        dsPan,dsMulti=None,None
        return []

  @staticmethod
  def writeMultispectralGeotiff( arrays,dsmulti,outname ):
    '''function writeMultispectralGeotiff( arrays,dsmulti,outname ):
    This function writes out a multispectral Geotiff (3 or 4 bands) 
    containing the Red,Green,Blue, and optionally, NIR bands. This 
    data should be passed into this function as a list[] of NumPy 
    arrays with 3 or 4 arrays. This function also takes in a GDAL 
    dataset object (dsmulti), which should have been created during
    the running of this command-line program as a 3 or 4-band 
    Geotiff image tile (a tile that is part of the larger resampled 
    multispectral Geotiff image file, resampled to match same 
    dimensions as panchromatic image file.)

    Args: 
      arrays (list): Python list[] holding Red,Green,Blue,NIR bands.
      dsmulti (osgeo.gdal.Dataset): GDAL dataset object for multispectral Geotiff.
      outname (str): output Geotiff filename.
    '''

    # if the multispectal geotiff already exists, delete it.
    if os.path.isfile(outname):
      os.remove(outname)
	  
    # write the multispectral geotiff
    nrows,ncols = arrays[0].shape
    driv = gdal.GetDriverByName('GTiff')
    dst = driv.Create( outname, ncols, nrows, len(arrays), gdal.GDT_Float32)
    dst.SetGeoTransform(dsmulti.GetGeoTransform())
    dst.SetProjection(dsmulti.GetProjection())
	
    # ------------------------------------------------------------------
    # in this method, we expect python list[] of 3 or 4 NumPy arrays,
    # in the following order: 
    #   (1) Red
    #   (2) Green
    #   (3) Blue
    #   (4) NIR (optional)
    # and we write the multispectral (multiband) geotiff with bands 
    # in the same order.
    # ------------------------------------------------------------------
	
    dst.GetRasterBand(1).WriteArray( arrays[0] )
    dst.GetRasterBand(2).WriteArray( arrays[1] )
    dst.GetRasterBand(3).WriteArray( arrays[2] )
	
    if len( arrays )>3:
      dst.GetRasterBand(4).WriteArray( arrays[3] )

    dst=None
    del dst

  @staticmethod
  def createMosaic( outname, outputTileNames ):
    '''function createMosaic( outname, outputTileNames ): 
    This static (class) method calls the gdal_merge.py Python script (which should be
    installed onto your system) to stitch-together a list of Geotiff image
    files into a single Geotiff. Hence, this function takes the list of Geotiff image
    files in outputTileNames, and passes them to gdal_merge.py such that they may
    be merged into a single Geotiff by a mosaic process. 

    Args:
      outname (str): Output string name of Geotiff mosaic.
      outputTileNames (list): List of strings for Geotiffs that should be mosaic together.
    '''
    mosaicCmd = gdalMergeLocation+' -o '+outname+' -of GTiff '
    mosaicCmdList = []
    mosaicCmdList.extend([mosaicCmd])
    mosaicCmdList.extend(outputTileNames)
    mosaicCmd = ' '.join(mosaicCmdList)
    proc = subprocess.Popen([mosaicCmd],shell=True,stdout=subprocess.PIPE)
    proc.wait()

  @staticmethod
  def resample( srcImageFilename, sourceDataset, dstDataset, outname, interp ): 
    '''function resample( srcImageFilename,sourceDataset,dstDataset,outname,interp):
    This function resamples a low-resolution multispectral Geotiff to larger 
    dimensions by means of the "interp" method (i.e. gdalconst.GRA_Cubic) passed
    into this function. Ultimately, this function resamples or resizes the  
    multispectral geotiff dataset, referred to by "sourceDataset" (and srcImageFilename)
    to the same dimensions as the panchromatic Geotiff image Geotiff file by 
    means of an interpolation method (i.e. bicubic resampling). 

    Args:
      srcImageFilename (str): source (low-res.) multispectral Geotiff filename. 
      sourceDataset (osgeo.gdal.Dataset): input multispectral GDAL dataset object.
      dstDataset (osgeo.gdal.Dataset): destination (high-res.) panchromatic dataset object. 
      outname (str): name of outputted resampled Geotiff
      interp (int): GDAL interpolation method (i.e. gdalconst.GRA_Cubic) 
    '''
      
    # get the "source" (i.e. low-res. multispectral) projection and geotransform
    srcProjection   = sourceDataset.GetProjection()
    srcGeotransform = sourceDataset.GetGeoTransform()
    srcNumRasters   = sourceDataset.RasterCount
    dstProjection   = dstDataset.GetProjection()
    dstGeotransform = dstDataset.GetGeoTransform()
    nrows           = dstDataset.RasterYSize
    ncols           = dstDataset.RasterXSize
    dst_fn          = outname
 
    # if the resampled-multispectral (3 or 4 band) Geotiff image file exists, delete it.
    if not os.path.isfile(outname):
      dst_ds = gdal.GetDriverByName('GTiff').Create( dst_fn, ncols,nrows,srcNumRasters,gdalconst.GDT_Float32) 
      dst_ds.SetGeoTransform(dstGeotransform)
      dst_ds.SetProjection(dstProjection)
      gdal.ReprojectImage(sourceDataset, dst_ds, srcProjection, dstProjection,interp)
      dst_ds=None
      del dst_ds
    return dst_fn

def getPath(name):
  '''function getPath(name):
  This function returns the path of an executable for command-line 
  tool on a Linux (UNIX-like) operating system by means of the 
  "which" command. For example, if the user passed-in the string
  "gdal_merge.py", and the full path on the system is 
  "/usr/bin/gdal_merge.py" , then this full-path string would 
  be returned.

  Args:
    name (str): name of command (i.e. gdal_merge.py) 
  Returns: 
    list: Output from "which" command for path of input name string.
  '''
  command = 'which '+name
  proc = subprocess.Popen([command],shell=True,stdout=subprocess.PIPE)
  proc.wait()
  return proc.stdout.readlines()

def usage():
  '''function usage(): 
  This function prints a usage message to the terminal. It can be invoked
  by running this script on the Linux (UNIX) command line with the -h or --help
  flags: 
    $ python pansharpen.py --help
  '''

  print '''
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
  '''
  sys.exit(1)

def version():
  print '''
   
    pansharpen.py: Version 1.0.0
    February 2019

    @Gerasimos Michalitsianos
    gerasimosmichalitsianos@gmail.com

  '''
  sys.exit(1)

def main():

  global gdalRetileLocation
  global gdalMergeLocation

  # ------------------------------------------------
  # declare empty strings for filename of:
  # (1) 1-band panchromatic Geotiff image file
  # (2) 3 or 4 band multispectral Geotiff image file
  # ------------------------------------------------

  multispectralGeotiff , panchromaticGeotiff = '' , '' 
  gdalMergeLocation , gdalRetileLocation = '' , ''

  try: 
    options,arguments = getopt.getopt(sys.argv[1:],'h:vpmsr:',['help','version',\
            'panchromatic=','multispectral=','gdalmerge=','gdalretile='])
  except getopt.GetoptError:
    usage()

  for option,argument in options:
    if option in ('-h','--h','--help'):
      usage()
    elif option in ('-p','--panchromatic'):
      panchromaticGeotiff = argument
    elif option in ('-v','--v','--version'):
      version()
    elif option in ('-m','--multispectral'):
      multispectralGeotiff = argument
    elif option in ('-s','--gdalmerge'):
      gdalMergeLocation = argument
    elif option in ('-r','--gdalretile'):
      gdalRetileLocation = argument
    else: pass

  # ---------------------------------------------------------------------
  # for both the gdal_merge.py,gdal_retile.py command-line programs that 
  # usually ship with GDAL, check if user passed them in at command line.
  # if they did not (i.e. its still an empty string), then attempt to 
  # get the path for each using "which" command ... if that fails, exit
  # the program.
  # ---------------------------------------------------------------------

  if gdalMergeLocation == '':  # gdal_merge.py location NOT passed in as command-line arg.
    gdalMergeLocation  = getPath('gdal_merge.py')
    if len(gdalMergeLocation)<1:
      print ' \n  gdal_merge.py command-line tool not passed-in at command line, and not found on system. Exiting ... \n '
      sys.exit(1)
    gdalMergeLocation  = gdalMergeLocation[0].strip()

  if gdalRetileLocation == '': # gdal_retile.py location NOT passed in as command-line arg. 
    gdalRetileLocation = getPath('gdal_retile.py')
    if len(gdalRetileLocation)<1:
      print ' \n  gdal_retile.py command-line tool not passed-in at command-line, and not found on system. Exiting ... \n' 
      sys.exit(1)
    gdalRetileLocation = gdalRetileLocation[0].strip()

  # ----------------------------------------------------------------------
  # make sure gdal_retile.py,gdal_merge.py paths have appropriate endings.
  # ----------------------------------------------------------------------

  if not gdalMergeLocation.endswith('gdal_merge.py') or not os.path.isfile(gdalMergeLocation): 
    print ' \n  not a valid path for gdal_merge.py: ', gdalMergeLocation , ' . Exiting ... \n '
    sys.exit(1)

  if not gdalRetileLocation.endswith('gdal_retile.py') or not os.path.isfile(gdalRetileLocation): 
    print ' \n  not a valid path for gdal_retile.py: ', gdalRetileLocation , ' . Exiting ... \n '
    sys.exit(1)
  
  print gdalRetileLocation
  print gdalMergeLocation
  sys.exit(1)

  # -------------------------------------------------
  # make sure both the panchromatic and multispectral 
  # Geotiff image filenames were passed-in via the 
  # the command-line.
  # -------------------------------------------------

  if panchromaticGeotiff == '' or multispectralGeotiff == '':
    print '  \n    Please pass in the names of a panchromatic (1-band) and multispectral (3 or 4 band) geotiff image files.'
    print '    This is done using the --panchromatic (or -p) and --multispectal (or -m) command-line flags.'
    usage()

  # -------------------------------------------------------
  # make sure both file do indeed exist on the file-system.
  # -------------------------------------------------------

  if not os.path.isfile( panchromaticGeotiff ): 
    print '  \n    Panchromatic geotiff image file does not exist: ' + panchromaticGeotiff + '. Exiting ... '
    usage()
  elif not os.path.isfile( multispectralGeotiff ): 
    print '  \n    Multispectral geotiff image file does not exist: ' + multispectralGeotiff + '. Exiting ... '
    usage()
  else: pass

  # ----------------------------------------------------
  # open up Panchromatic and Multispectral Geotiff image 
  # files as GDAL objects in Python memory. 
  # ----------------------------------------------------

  dsPan   = gdal.Open( panchromaticGeotiff , gdal.GA_ReadOnly )
  dsMulti = gdal.Open( multispectralGeotiff, gdal.GA_ReadOnly )

  # ----------------------------------------------------
  # make sure Panchromatic Geotiff has only 1 band, and 
  # Multispectral Geotiff image file has 3 or 4 bands.
  # ----------------------------------------------------

  if dsPan.RasterCount != 1:
    print '  \n    Panchromatic Geotiff image file: ' + panchromaticGeotiff + ' should have ONE single band. Exiting ... '
    sys.exit(1) 

  if (dsMulti.RasterCount != 3) and ( dsMulti.RasterCount != 4 ):
    print '  \n    Multispectral Geotiff image file: ' + multispectralGeotiff + ' should have 3 or 4 bands. Exiting ... '
    sys.exit(1) 

  # -----------------------------------------------------
  # resample 3 or 4 band Multispectral Geotiff image file
  # using bicubic resampling.
  # -----------------------------------------------------

  if multispectralGeotiff.endswith('.TIF'):
    resampledMultispectralGeotiffFilename = multispectralGeotiff.replace('.TIF', '_RESAMPLED.TIF')
  elif multispectralGeotiff.endswith('.tif'): 
    resampledMultispectralGeotiffFilename = multispectralGeotiff.replace('.tif', '_RESAMPLED.TIF')
  else: 
    print '  \n    Multispectral Geotiff image file: ' + multispectralGeotiff + ' should have .TIF or .tif extension. Exiting ... '
    sys.exit(1)

  # --------------------------------------------------
  # perform bicubic multispectral resampling of 3 or 4
  # band Geotiff image file.
  # --------------------------------------------------

  resampledMultispectralGeotiffFilename = pansharpen.resample(
    multispectralGeotiff, 
    dsMulti, 
    dsPan, 
    resampledMultispectralGeotiffFilename, 
    gdalconst.GRA_Cubic
  )

  # ---------------------------------------------------
  # At this point, exit the program if we do not have 
  # the resampled 4-band Geotiff image file.
  # ---------------------------------------------------

  if not os.path.isfile(resampledMultispectralGeotiffFilename): 
    print '  \n    Multispectral Geotiff image file: ' + multispectralGeotiff + ' . FAILURE to resample. Exiting .... '
    sys.exit(1)

  # ---------------------------------------------------------------------------
  # create output Geotiff filenames for 4 different multispectral (3 or 4 band)
  # Geotiff image filenames: 
  #  (1) pan-sharpened FIHS (Fast Intensity Hue Saturation) Geotiff file.
  #  (2) pan-sharpened Brovey Geotiff file.
  #  (3) PCA pan-sharpened Geotiff file. 
  #  (4) Wavelet pan-sharpened Geotiff file.
  # ---------------------------------------------------------------------------

  outnameFIHS    = resampledMultispectralGeotiffFilename.replace(
    '_RESAMPLED.TIF', '_panSharpenedFIHS.tif')
  outnameBROVEY  = resampledMultispectralGeotiffFilename.replace(
    '_RESAMPLED.TIF', '_panSharpenedBROVEY.tif')
  outnamePCA = resampledMultispectralGeotiffFilename.replace(
    '_RESAMPLED.TIF', '_panSharpenedPCA.tif')
  outnameWAVELET = resampledMultispectralGeotiffFilename.replace(
    '_RESAMPLED.TIF', '_panSharpenedWAVELET.tif')
  
  # -----------------------------------------------------------------
  # if any of the pan-sharpened Brovey,FIHS, or PCA Geotiff files 
  # already exist, delete them.
  # -----------------------------------------------------------------

  if os.path.isfile( outnameFIHS ): 
    os.remove( outnameFIHS )
  elif os.path.isfile( outnameBROVEY ):
    os.remove( outnameBROVEY )
  elif os.path.isfile( outnamePCA ):
    os.remove( outnamePCA )
  elif os.path.isfile( outnameWAVELET ):
    os.remove( outnameWavelet )
  else: pass 

  # ---------------------------------------------------------------
  # tile up resampled (very large) multispectral Geotiff image file
  # ---------------------------------------------------------------

  # define number of tiles in one particular direction
  dim=5.0
  nrows,ncols = dsPan.RasterYSize, dsPan.RasterXSize
  xTileDim,yTileDim = int(ncols/dim) , int(nrows/dim)

  # create tile command and then subsequently run it. tiles go to same directory as Pan,MS Geotiffs
  pythonInterpreter  = sys.executable 
  targetDirectory = os.path.dirname( panchromaticGeotiff )
  multispectralTilesListCSV = os.path.join( 
    targetDirectory, 'multispectral_Tilenames.CSV')
  panchromaticTilesListCSV  = os.path.join( 
    targetDirectory, 'panchromatic_Tilenames.CSV' )

  # create commands in which to tile up both the Panchromatic and Multispectral Geotiff image files.
  tileCmdMS = pythonInterpreter + ' ' + gdalRetileLocation + ' -ps ' + str(xTileDim) + ' ' + \
    str(yTileDim) + ' -targetDir ' + targetDirectory + ' -csv ' + \
    os.path.basename(multispectralTilesListCSV) + ' ' + resampledMultispectralGeotiffFilename
  tileCmdPan = pythonInterpreter + ' ' + gdalRetileLocation + ' -ps ' + str(xTileDim) + ' ' + \
    str(yTileDim) + ' -targetDir ' + targetDirectory + ' -csv ' + \
    os.path.basename(panchromaticTilesListCSV) + ' ' + panchromaticGeotiff

  # produce the tiles for the Panchromatic and [Resampled] Multispectral Geotiff image files.
  procTileMS  = subprocess.Popen([tileCmdMS],shell=True,stdout=subprocess.PIPE)
  procTileMS.wait()

  procTilePan = subprocess.Popen([tileCmdPan],shell=True,stdout=subprocess.PIPE)
  procTilePan.wait()

  # make sure tiles were successfully created ... this would be evident by the CSV with Geotiff tile filenames
  if not os.path.isfile( multispectralTilesListCSV ): 
    print '  \n    Unable to produce tiles for: ' + resampledMultispectralGeotiffFilename + ' . FAILURE to tile. Exiting .... '
    sys.exit(1)

  if not os.path.isfile( panchromaticTilesListCSV ): 
    print '  \n    Unable to produce tiles for: ' + panchromaticGeotiff + ' . FAILURE to tile. Exiting .... '
    sys.exit(1)
 
  # parse the CSVs to get the names of the 
  # Multispectral (resampled) and panchromatic Geotiff files

  linesTilesMS  = open( multispectralTilesListCSV , 'r').readlines()
  linesTilesPan = open( panchromaticTilesListCSV  , 'r').readlines()

  # --------------------------------------------------------
  # initialize Python lists to hold pan-sharpened tile-names 
  # --------------------------------------------------------

  outputTileNamesBrovey, outputTileNamesFIHS, outputTileNamesPCA, outputTileNamesWavelet, tilenamesMS, tilenamesPan \
          = ([] for i in range(6))

  # ----------------------------------
  # iterate through Geotiff tile-names
  # ----------------------------------

  for msTileNameLine,panTileNameLine in zip( linesTilesMS,linesTilesPan ): 
  
    # -------------------------------------------------------
    # get names of 3 or 4-band bicubic-resampled Geotiff tile
    # as well as the panchromatic Geotiff image tile filename
    # -------------------------------------------------------

    tilenameMS  = os.path.join( targetDirectory, msTileNameLine.split(';')[0] )
    tilenamePan = os.path.join( targetDirectory, panTileNameLine.split(';')[0] )

    # --------------------------------------------------------
    # create pan-sharpening object
    # --------------------------------------------------------

    imgFusion = pansharpen( tilenamePan, tilenameMS )
	
    # ---------------------------------------------------
    # produce lists[] of pan-sharpened arrays. Each list
    # contains 4 NumPy arrays: 
    #   (1) Red
    #   (2) Green
    #   (3) Blue
    #   (4) NIR (optional)
    # ---------------------------------------------------

    sharpenedFIHS    = imgFusion.pansharpenFIHS()
    sharpenedBrovey  = imgFusion.pansharpenBrovey()
    sharpenedPCA     = imgFusion.pansharpenPCA()
    sharpenedWavelet = imgFusion.pansharpenWavelet()

    # ---------------------------------------------------
    # establish Geotiff out-names for FIHS,Brovey,PCA
    # ---------------------------------------------------

    outnameTileBrovey  = tilenameMS.replace(
            '.tif','_PanSharpenedBrovey.tif')
    outnameTileFIHS    = tilenameMS.replace(
            '.tif','_PanSharpenedFIHS.tif')
    outnameTilePCA = tilenameMS.replace(
            '.tif','_PanSharpenedPCA.tif')
    outnameTileWavelet = tilenameMS.replace(
            '.tif','_PanSharpenedWavelet.tif')

    # ---------------------------------------------
    # write both pan-sharpened multispectral images 
    # to geotiffs (for one tile sector)
    # ---------------------------------------------
	
    dsmulti=gdal.Open( tilenameMS )
    pansharpen.writeMultispectralGeotiff( sharpenedFIHS, dsmulti, outnameTileFIHS )
    pansharpen.writeMultispectralGeotiff( sharpenedBrovey, dsmulti, outnameTileBrovey )
    pansharpen.writeMultispectralGeotiff( sharpenedPCA, dsmulti, outnameTilePCA )
    pansharpen.writeMultispectralGeotiff( sharpenedWavelet, dsmulti, outnameTileWavelet )
    dsmulti=None
    del dsmulti
	
    outputTileNamesBrovey.append(outnameTileBrovey)
    outputTileNamesFIHS.append( outnameTileFIHS)
    outputTileNamesPCA.append( outnameTilePCA )
    outputTileNamesWavelet.append( outnameTileWavelet)

    tilenamesMS.append(tilenameMS)
    tilenamesPan.append(tilenamePan)
   
  # ----------------------------------------------------------
  # perform mosaic process for 3 or 4 band Geotiff image files
  # ----------------------------------------------------------
  
  pansharpen.createMosaic( outnameBROVEY, outputTileNamesBrovey )
  pansharpen.createMosaic( outnameFIHS, outputTileNamesFIHS )
  pansharpen.createMosaic( outnamePCA, outputTileNamesPCA )
  pansharpen.createMosaic( outnameWAVELET, outputTileNamesWavelet )

  # clean up tiles 
  for tilePan,tileMS,tileBrovey,tileFIHS,tilePCA,tileWavelet in \
          zip(tilenamesPan,tilenamesMS,outputTileNamesBrovey,outputTileNamesFIHS,outputTileNamesPCA,outputTileNamesWavelet):
    os.remove(tileBrovey)
    os.remove(tileFIHS)
    os.remove(tileMS)
    os.remove(tilePan)
    os.remove(tilePCA)
    os.remove(tileWavelet)

  # clean up CSV files
  os.remove(multispectralTilesListCSV)
  os.remove(panchromaticTilesListCSV)
  os.remove(resampledMultispectralGeotiffFilename)

if __name__ == '__main__':
  main()
