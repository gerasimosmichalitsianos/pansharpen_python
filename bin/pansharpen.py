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

  def __init__(self,outputDirectory,red,green,blue,pan,NIR=None):
    '''function __init__(self,pan,multi):
    This is a constructor method to create an object of this class.
    Instantiate an object of this class with 4 or 5 input NumPy arrays
    holding data for the following bands of a satellite image: 
      (1) Red band 
      (2) Green band 
      (3) blue band 
      (4) panchromatic band (grey-scale)
      (5) NIR band , near-infrared band (optional)
    Returns: 
      pancharpen: an object of the pansharpen class.

    '''
    self.outDir = outputDirectory
    self.red = red
    self.green = green
    self.blue = blue
    self.pan = pan

    if not NIR is None: 
      self.NIR = NIR
    else: 
      self.NIR = None

  def pansharpenPCA(self): 
    '''function pansharpenPCA(self):
    This is an instance method that returns a Python list of 3 or 4
    NumPy arrays containing the pan-sharpened Red,Green,Blue, and 
    optionally, NIR bands. These bands will have been created using 
    the principal component analysis (PCA) pan-sharpening.
    
    Returns: 
      list: Python list[] containing 3 or 4 NumPy arrays using PCA method.
    '''
    nrows,ncols = self.pan.shape
    if self.NIR is not None:
      image = np.zeros( (nrows,ncols,4),dtype=np.float32)
      image[:,:,0] = self.red
      image[:,:,1] = self.green
      image[:,:,2] = self.blue
      image[:,:,3] = self.NIR
    else:
      image = np.zeros( (nrows,ncols,3) , dtype=np.float32)
      image[:,:,0] = self.red
      image[:,:,1] = self.green
      image[:,:,2] = self.blue

    m,n,d=image.shape
    nanlocs = np.where(self.pan==0.0)
    image   = np.where( image == 0.0, np.nan, image )
    image   = np.ma.array(image, mask=np.isnan(image))

    M = np.reshape( image, (m*n,d))
    [ PCAData , PCAMap ] = pansharpen.multiPCA( M,d )
    PCAData = np.reshape( PCAData , (m,n,d))
    F = PCAData

    P = ( self.pan - self.pan.mean() ) * (np.std(F) / self.pan.std() )  + F.mean()
    F[:,:,0] = P
    F = pansharpen.inversePCA( PCAMap['M'] , np.reshape( F, (m*n,d)), PCAMap['mean'])
    sharpened = np.reshape(F,(m,n,d))

    outred   = sharpened[ :,:,0 ]
    outgreen = sharpened[ :,:,1 ]
    outblue  = sharpened[ :,:,2 ]

    if self.NIR is not None:
      outNIR = sharpened[ :,:,3 ]
      return [ outred, outgreen, outblue, outNIR ]
    else:
      return [ outred, outgreen, outblue ]

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
    if self.NIR is None:

      L = ( self.red+self.green+self.blue ) / 3.0
      redsharp   = self.red   + ( self.pan   - L )
      greensharp = self.green + ( self.green - L )
      bluesharp  = self.blue  + ( self.blue  - L )
      return [redsharp,greensharp,bluesharp]

    elif self.NIR is not None: 

      L = ( self.red+self.green+self.blue+self.NIR ) / 4.0
      redsharp   = self.red   + ( self.pan   - L )
      greensharp = self.green + ( self.green - L )
      bluesharp  = self.blue  + ( self.blue  - L )
      NIRsharp   = self.NIR   + ( self.NIR   - L )
      return [redsharp,greensharp,bluesharp,NIRsharp]

    else: 
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

    # set number of output bands (3 or 4, depending on if NIR is passed-in)
    if self.NIR is not None:
      d = 4
    else: 
      d = 3

    # get 2D dimensions of output imagery
    nrows,ncols = self.pan.shape
    
    # create 3D data cube to perform wavelet pan-sharpening
    if d == 3: 
      image = np.zeros((nrows,ncols,d),dtype=np.float16)
      image[:,:,0] = self.red
      image[:,:,1] = self.green
      image[:,:,2] = self.blue
    elif d == 4:
      image = np.zeros((nrows,ncols,d),dtype=np.float16)
      image[:,:,0] = self.red
      image[:,:,1] = self.green
      image[:,:,2] = self.blue
      image[:,:,3] = self.NIR
    else: 
      return []

    level = 0
    wavelet_type = 'haar'

    coeffs = pywt.wavedec2( self.pan, wavelet=wavelet_type, level=level)
    panvec,coeff_slices,coeff_shapes = pywt.ravel_coeffs(coeffs)
    reconstvec = np.tile( panvec.T, (d,1)).T

    n=panvec.shape[0]
    lowresvec = np.zeros((n,d),dtype=np.float16)

    for band in range(d):
      lowresCoeffs = pywt.wavedec2( image[:,:,band], wavelet=wavelet_type, level=level)
      lowresArr,arrSlices = pywt.coeffs_to_array(lowresCoeffs)
      lowresvec[:,band] = np.reshape(lowresArr,(nrows*ncols,))

    for j in range( 0 , coeff_shapes[0][0] * coeff_shapes[0][1] ):
      reconstvec[ j,: ] = lowresvec[j,:]

    sharpened = np.zeros( (nrows,ncols,d), dtype=np.float16)
    for band in range(d):
      p = np.reshape( reconstvec[:,band], (nrows,ncols))
      fcoeffs = pywt.wavedec2(p,wavelet_type,level=level)
      out=pywt.waverec2(fcoeffs,wavelet_type)
      sharpened[:,:,band] = out

    redsharp   = sharpened[:,:,0]
    greensharp = sharpened[:,:,1]
    bluesharp  = sharpened[:,:,2]
   
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
    with warn.catch_warnings():
      warn.filterwarnings('ignore',category=RuntimeWarning)

      if self.NIR is None: 

        redsharp   = np.multiply( np.true_divide( 
          self.red, self.red+self.green+self.blue)   , self.pan )
        greensharp = np.multiply( 
          np.true_divide( self.green, self.red+self.green+self.blue) , self.pan )
        bluesharp  = np.multiply( 
          np.true_divide( self.blue, self.red+self.green+self.blue)  , self.pan )
        return [redsharp,greensharp,bluesharp]

      elif self.NIR is not None:
        NIR   = self.NIR
        redsharp   = np.multiply( np.true_divide( 
          self.red, self.red+self.green+self.blue+self.NIR)   , self.pan )
        greensharp = np.multiply( np.true_divide( 
          self.green, self.red+self.green+self.blue+self.NIR) , self.pan )
        bluesharp  = np.multiply( np.true_divide( 
          self.blue, self.red+self.green+self.blue+self.NIR)  , self.pan )
        NIRsharp   = np.multiply( np.true_divide( 
          self.NIR, self.red+self.green+self.blue+self.NIR)   , self.pan )
        return [redsharp,greensharp,bluesharp,NIRsharp]
      else:
        return []

  @staticmethod
  def writeMultispectralGeotiff( arrays,dsmulti,outname ):
    '''function writeMultispectralGeotiff( arrays,dsmulti,outname ):
    This function writes out a multispectral Geotiff (3 or 4 bands) 
    containing the Red,Green,Blue, and optionally, NIR bands. This 
    data should be passed into this function as a list[] of NumPy 
    arrays with 3 or 4 arrays. This function also takes in a GDAL 
    dataset object (dsmulti), which should have been created during
    the running of this command-line program. 

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

def usage():
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
  '''
  sys.exit(1)

def main():

  # ------------------------------------------------
  # declare empty strings for filename of:
  # (1) 1-band panchromatic Geotiff image file
  # (2) 3 or 4 band multispectral Geotiff image file
  # ------------------------------------------------

  multispectralGeotiff , panchromaticGeotiff = '' , '' 

  try: 
    options,arguments = getopt.getopt(sys.argv[1:],'h:p:m:',['help','version',\
      'panchromatic=','multispectral=',])
  except getopt.GetoptError:
    usage()

  for option,argument in options:
    if option in ('-h','--h','--help'):
      usage()
    elif option in ('-p','--panchromatic'):
      panchromaticGeotiff = argument
    elif option in ('-m','--multispectral'):
      multispectralGeotiff = argument
    else: pass

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
 
  # ----------------------------------
  # clear up datasets after resampling
  # ----------------------------------

  dsMulti,dsPan = None,None
  del dsMulti,dsPan

  # -----------------------------------------------------
  # At this point, exit the program if we do not have 
  # the resampled 4-band (or 3-band) Geotiff image file.
  # -----------------------------------------------------

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

  # -------------------------------------------------------------
  # use Python/NumPy's memorymap function(s) to store our arrays
  # for Red,Green,Blue,NIR(optional), and Panchromatic bands 
  # as binary files on disk.
  # -------------------------------------------------------------

  # create GDAL dataset objects for panchromatic geotiff and resampled 
  # multispectral (3 or 4 band) geotiff
  dsPan   = gdal.Open( panchromaticGeotiff )
  dsMulti = gdal.Open( resampledMultispectralGeotiffFilename )

  # create output directory to hold .dat files (binary files)
  outputDir = os.path.dirname( panchromaticGeotiff)
  nrows,ncols = dsPan.RasterYSize, dsPan.RasterXSize

  red   = dsMulti.GetRasterBand(1).ReadAsArray().astype(float)
  green = dsMulti.GetRasterBand(2).ReadAsArray().astype(float)
  blue  = dsMulti.GetRasterBand(3).ReadAsArray().astype(float)
  pan   = dsPan.GetRasterBand(1).ReadAsArray().astype(float)

  if dsMulti.RasterCount>3:
    NIR   = dsMulti.GetRasterBand(4).ReadAsArray().astype(float)
    imgFusion  = pansharpen( outputDir, red,green,blue,pan,NIR)
  else: 
    imgFusion  = pansharpen( outputDir, red,green,blue,pan )

  sharpenedFIHS    = imgFusion.pansharpenFIHS()
  pansharpen.writeMultispectralGeotiff( sharpenedFIHS,dsMulti,outnameFIHS )
  del sharpenedFIHS

  sharpenedBrovey  = imgFusion.pansharpenBrovey()
  pansharpen.writeMultispectralGeotiff( sharpenedBrovey,dsMulti,outnameBROVEY )
  del sharpenedBrovey

  sharpenedWavelet = imgFusion.pansharpenWavelet()
  pansharpen.writeMultispectralGeotiff( sharpenedWavelet,dsMulti,outnameWAVELET )
  del sharpenedWavelet 

  sharpenedPCA     = imgFusion.pansharpenPCA()
  pansharpen.writeMultispectralGeotiff( sharpenedPCA,dsMulti,outnamePCA )
  del sharpenedPCA
  
  # clean-up resampled MS geotiff ... we don't need it anymore
  # ----------------------------------------------------------

  os.remove( resampledMultispectralGeotiffFilename )

if __name__ == '__main__':
  main()
