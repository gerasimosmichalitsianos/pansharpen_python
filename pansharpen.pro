pro pansharpenimage

  ;
  ; This IDL script is an example of how to use 
  ; ENVI functions to read Geotiff files in, and
  ; perform a Gram-Schmidt pan-sharpening. In this
  ; process, the multispectral Geotiff file 
  ; is pan-sharpened to the same dimensions as the 
  ; input panchromatic Geotiff image. To use this 
  ; script, both IDL and ENVI must be installed. 
  ; 
  ; usage
  ; $ idl -e "pansharpen" -args panfname.tif multifname.tif 
  ; 
  ; @author: 
  ; Gerasimos Michalitsianos 
  ; NASA/GSFC, Science Systems and Applications, Inc. 
  ; May 2015 
  ;

  ; go to appropriate directory location for files, start ENVI
  envi 
  
  ; get our Geotiff filenames 
  ; (1) panchromatic geotiff file (1 band)
  ; (2) multispectral geotiff file (4 bands, RGB, NIR)
  
  args = command_line_args() 
  usageMessage = "$ idl -e 'pansharpen' panfname.tif multifname.tif"  

  panfname = args[0]
  multifname = args[1]
  
  if panfname.endswith('.tif') ne 1 then begin 
    print, 'Not a Geotiff file: ' + panfname 
    return 
  endif 
  
  if multifname.endswith('.tif') ne 1 then begin 
    print, 'Not a Geotiff file: ' + multifname 
    return
  endif

  ; open up the coarse multispectral file, get its image and map parameters
  envi_open_data_file, msfname, r_fid = fidlow, /tiff
  if (fidlow eq -1) then return 
  envi_file_query, fidlow, nb = nblow, dims = dimslow

  ; open up the high-resolution panchromatic image file, get its image and map parameters
  envi_open_data_file, panfname, r_fid = fidhigh, /tiff
  if (fidhigh eq -1) then return 
  envi_file_query, fidhigh, nb = nbhigh, dims = dimshigh

  ; create an outfile name
  outname = msfname.replace('.tif','_PanSharpenedGramSchmidt.tif')
  outname = outname.trim()
  pos=[0,1,2,3]

  ; now we are going to perform some pan-sharpening
  envi_doit, 'envi_gs_sharpen_doit', $

    fid = fidlow, $
    dims = dimslow, $
    pos = pos, $
    method=0, $
    out_name = outname, $
    interp=2, $
    hires_fid = fidhigh, $
    hires_pos = [0], $ 
    hires_dims = dimshigh

end 
