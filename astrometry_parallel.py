import numpy as np
import pyfits as pf
import os
import glob
import time
from multiprocessing import Pool
import sys, getopt

###############################################################################
#modify these and only these variables:

period='20130626' #period for which the list will be generated

ncores=28 #number of cores used to combine the darks. A high fraction of the total is not recomended, since some process use multiple cores.

main_path='/paula3/QUEST/data/'  #where are the daily image folders located

defa='/paula3/QUEST/scripts/default_files'

mask_path='/paula3/QUEST/data/Masks/'

###############################################################################

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"p:n:d:defa:")
###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-p':
        period=a
    elif o == '-n':
        ncores=a
    elif o == '-d':
        main_path=a
    elif o == '-defa':
        defa=a


period=str(period) #period for which the list will be generated

ncores=int(ncores) #number of cores used in the reduction process

main_path=str(main_path)


###############################################################################


start = time.clock()


list_gain=[5.405,5.128,5.988,5.154,5.0,8.474,5.263,6.211,4.524,5.0,5.376,0.0,5.747,5.102,0.0,4.926,5.347,5.714,0.0,4.926,5.347,5.524,5.154,5.586,4.608,15.625,5.524,5.0,0.0,5.291,4.739,5.617,6.134,5.434,0.0,5.555,5.235,5.952,5.988,5.128,4.878,5.780,5.882,5.649,5.555,6.134,5.347,6.756,5.405,4.950,5.347,5.076,6.369,5.917,4.739,5.617,5.154,6.410,6.060,5.617,5.050,5.0,0.0,4.444,5.747,5.128,5.376,5.681,5.102,5.347,5.747,5.0,5.813,5.649,5.917,5.714,8.547,6.993,6.993,5.555,5.988,8.849,4.830,5.747,6.024,6.622,6.289,6.134,4.807,9.345,0.0,7.692,5.617,5.494,5.847,5.988,5.376,5.050,5.988,5.524,4.807,4.950,5.649,5.263,4.672,5.524,5.154,5.181,5.780,5.813,7.751,5.405]
list_gain=np.array(list_gain)


chips_name=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']

chips_name=np.array(chips_name)


img_type='s.'

list_folders=sorted(glob.glob(main_path+period+'*/*'+img_type+'*/reduced'))
#list_folders=np.core.defchararray.replace(list_folders,main_path,'')
#list_folders=np.array(list_folders)


print list_folders

seeing=3.5

year=period[0:4]

def run_astrometry(sci_folder):

    os.system('rm '+sci_folder+'/R*.new.fits >& /dev/null')
    os.system('rm '+sci_folder+'/R*_sel.cat >& /dev/null')
    os.system('rm '+sci_folder+'/R*cat.fits >& /dev/null')
    os.system('rm '+sci_folder+'/R*cat >& /dev/null')
    os.system('rm '+sci_folder+'/R*phot.fits >& /dev/null')

    rlist=sorted(glob.glob(sci_folder+'/R_*fits'))
    #print rlist

    if len(rlist)>0:

	list_chip=[x[-8:-5] for x in rlist]

    	rlist=np.core.defchararray.replace(rlist,'.fits','')
    	list_chip=np.array(list_chip)
	#print rlist

        if (os.path.exists(sci_folder+'/tmp')==False):
            os.system('mkdir '+sci_folder+'/tmp')
            print "new tmp folder created in %s " % (sci_folder)


	for i, item in enumerate(rlist):
            #name of the input image
            input_img=rlist[i]+'.fits'
            #print input_img

            sel_input=rlist[i]+'_sel.cat'

            #name of the output image with the astrometry solved
            wcs=rlist[i]+'.new.fits'

            #name of the output catalog to be used in the calibration
            catalog=rlist[i]+'.cat.fits'

	    #print input_img,sel_input,wcs,catalog

            nchip=np.where(list_chip[i]==chips_name)
            gain=list_gain[nchip][0]

	    print input_img,sel_input,wcs,catalog,gain

            if os.path.isfile(wcs): os.system('rm '+wcs)
            if os.path.isfile(catalog): os.system('rm '+catalog)
            if os.path.isfile(sel_input): os.system('rm '+sel_input)

            #we check if the image is good enough to calculate the astrometry

            #we run sextractor to see whether there are real deteccions in the image
            cm='sex '+input_img+' -c '+defa+'/quest.sex -CATALOG_NAME '+sel_input+' -CATALOG_TYPE FITS_1.0 -PARAMETERS_NAME '+defa+'/selec_chips.param -WEIGHT_TYPE NONE -GAIN '+str(gain)+' -SEEING_FWHM '+str(seeing)+' -FILTER_NAME '+defa+'/gauss_3.0_7x7.conv -STARNNW_NAME '+defa+'/default.nnw -PHOT_APERTURES 7.0'
            os.system(cm)

            #we check the sextractor results
            arch=pf.open(sel_input)
            datos=arch[1].data
            flux=datos['FLUX_APER']


            if len(flux)>0:

                flag=datos['FLAGS']
                frad=datos['FLUX_RADIUS']
                backg=datos['BACKGROUND']
                fwhm=datos['FWHM_IMAGE']
                clas=datos['CLASS_STAR']
                n=np.where((flag==0) & (fwhm>=2.5) & (fwhm<=7.0))
                flux=flux[n]
                backg=backg[n]
                fwhm=fwhm[n]
                clas=clas[n]

                if len(clas)>0:

                    diff=flux-backg
                    m=np.where((diff>10000.0) &(clas>0.5))
                    fwhm=fwhm[m]
                    l=len(diff)

                    #we select images with a considerable amount of sources detected
                    if l>19:
                        ss=np.median(fwhm)
                        sg=ss*0.882
                        seeing_chip=sg

                        #we select images with good seeing
                        if seeing_chip>=1.0:
                            #we run the astrometry
                            aa=pf.open(input_img)
                            head=aa[0].header
                            CHIPRA=head['CHIP-RA']
                            CHIPDEC=head['CHIP-DEC']

                            #we run astrometry.net using a gess for the ra, dec and pixel sclae
                            os.system(' /usr/local/astrometry/bin/solve-field  '+input_img+'  --temp-dir '+sci_folder+'/tmp  --overwrite --no-plots  --scale-units arcsecperpix --scale-low 0.88 --scale-high 0.89 --no-verify --cpulimit 120 --radius 5 --ra '+CHIPRA+' --dec '+CHIPDEC)

                            #if the astrometry is succesful we change the name of the output and run sextractor again to generate the catalogs to be used in the calibration
                            if os.path.isfile(rlist[i]+'.new'):

                                os.system('mv '+rlist[i]+'.new '+wcs)

                                os.system('rm '+rlist[i]+'*.png >& /dev/null')
                                os.system('rm '+rlist[i]+'*.xyls >& /dev/null')
                                os.system('rm '+rlist[i]+'*.rdls >& /dev/null')
                                os.system('rm '+rlist[i]+'*.axy >& /dev/null')
                                os.system('rm '+rlist[i]+'*.solved >& /dev/null')
                                os.system('rm '+rlist[i]+'*.match >& /dev/null')
                                os.system('rm '+rlist[i]+'*.kmz >& /dev/null')
                                os.system('rm '+rlist[i]+'*.corr >& /dev/null')
				os.system('rm '+rlist[i]+'*.wcs >& /dev/null')

                                #we run sextractor
	
				var_map=mask_path+year+'/'+'VarMap_'+year+'_'+list_chip[i]+'.fits'
    				flag_img=mask_path+year+'/'+'Mask_'+year+'_'+list_chip[i]+'.fits'
                                #print flag_img
				if os.path.isfile(var_map):
				    cmd='sex '+wcs+' -c '+defa+'/apphot.sex -CATALOG_TYPE FITS_1.0 -CATALOG_NAME '+catalog+' -PARAMETERS_NAME '+defa+'/apphot.param -FLAG_TYPE AND  -FLAG_IMAGE  '+flag_img+'   -WEIGHT_TYPE MAP_VAR  -WEIGHT_IMAGE '+var_map+' -GAIN '+str(gain)+' -SEEING_FWHM '+str(seeing_chip)+' -FILTER_NAME '+defa+'/gauss_3.0_7x7.conv -STARNNW_NAME '+defa+'/default.nnw'

                                    os.system(cmd)
				    #print cmd
				    print "running sextractor"
				os.system('rm '+sci_folder+'/tmp/*')


            os.system('rm '+sel_input)


#we run the astrometry using miltiprocessing
<pool = Pool(ncores) #ncores cores are used

results = pool.map(run_astrometry,list_folders)

pool.close()
pool.join()



elapsed = (time.clock() - start)
print elapsed
