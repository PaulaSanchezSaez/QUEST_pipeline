#code to reduce the science images
import time
import numpy as np
import pyfits as pf
import os
import glob
from astropy import units as u
import ccdproc as ccd
from multiprocessing import Pool
from skimage.restoration import inpaint
import sys, getopt

###############################################################################
#modify these and only these variables:

period='201511' #period for which the list will be generated

ncores=28 #number of cores used in the reduction process

mid_epoch=20151115 #epoch that defines the end of the 1st period and the start of the second period (for the flat correction)

main_path='/paula3/QUEST/data/'  #where are the daily image folders located

###############################################################################
# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"p:n:m:d:")
###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-p':
        period=a
    elif o == '-n':
        ncores=a
    elif o == '-m':
        mid_epoch=a
    elif o == '-d':
        main_path=a




period=str(period) #period for which the list will be generated

ncores=int(ncores) #number of cores used in the reduction process

mid_epoch=int(mid_epoch) #epoch that defines the end of the 1st period and the start of the second period (for the flat correction)

main_path=str(main_path)



###############################################################################


start = time.clock()

img_type='s.'

list_sci=sorted(glob.glob(main_path+period+'*/*'+img_type+'*/*'+img_type+'*fits')) #list with the fits images available

list_sci=np.core.defchararray.replace(list_sci,main_path,'')

list_sci=np.array(list_sci)

list_epoch=[x[0:8] for x in list_sci] #list with the epoch for every image

list_epoch=np.array(list_epoch)

list_folder=[x[0:26] for x in list_sci] #list with the folder for every image

list_folder=np.array(list_folder)

list_name=[x[27:51] for x in list_sci] #list with the name for every image

list_name=np.array(list_name)

list_chips=[x[43:46] for x in list_sci] #list with the name for every image

list_chips=np.array(list_chips)


epochs=sorted(glob.glob(main_path+period+'*')) #array with the available epochs

epochs=np.core.defchararray.replace(epochs,main_path,'')

epochs=np.array(epochs)


folders=sorted(glob.glob(main_path+period+'*/*s.*')) #array with the available science image folders

folders=np.core.defchararray.replace(folders,main_path,'')

folders=np.array(folders)

month=period[0:6]


print "epochs for which the masterdark will be generated = ", epochs


chips_name=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']

chips_name=np.array(chips_name)



def red_img(sci_img,sci_folder,chip,day):
#function to substract de dark from the flat

    output_name='R_'+sci_img

    print output_name

    #getting info from the science image
    a=pf.open(main_path+sci_folder+'/'+sci_img)
    head=a[0].header
    ext=head['EXPTIME']

    a.close()
    del a

    #defining the exposure time
    if (ext>56.0) and (ext<64.0):  exptime='60s'
    elif (ext>176.0) and (ext<184.0):  exptime='180s'

    #defining the flat period
    if int(day)<=mid_epoch: flat_period=month+'1st'
    elif int(day)>mid_epoch: flat_period=month+'2nd'

    #defining the dark for the science image
    dark_name=main_path+day+'/darks/'+'Dark_'+exptime+'_'+chip+'.fits'

    #defining the flat for the science image
    flat_name=main_path+'MasterFlats/'+flat_period+'/'+'Flat_'+flat_period+'_'+chip+'.fits'

    #defining the mask for the science image
    #mask_name=main_path+'Masks/Mask_old/'+'Mask_'+chip+'.fits'

    if (os.path.exists(main_path+sci_folder+'/reduced/')==False):
        os.system('mkdir '+main_path+sci_folder+'/reduced/')
        print "new reduced folder created in %s " % (sci_folder)


    if os.path.isfile(main_path+sci_folder+'/reduced/'+output_name):
        os.system('rm  '+main_path+sci_folder+'/reduced/'+output_name)
        print "image deleted %s" % (main_path+sci_folder+'/reduced/'+output_name)


    if os.path.isfile(dark_name) and os.path.isfile(flat_name):
    	#the image is trimmed
    	timg=ccd.trim_image(ccd.CCDData.read(main_path+sci_folder+'/'+sci_img,unit="adu"),fits_section='[5:595,5:2395]')
    	print "$$$$$$$$$$  trim done $$$$$$$$$$$$$$$", sci_img

    	#the darks is substracted
    	sdimg=ccd.subtract_dark(timg,ccd.CCDData.read(dark_name,unit="adu"),exposure_time='EXPTIME',exposure_unit=u.second,scale=False)
    	print "$$$$$$$$$$  subs dark done $$$$$$$$$$$$$$$", sci_img

    	#the flat is applied
    	fcimg=ccd.flat_correct(sdimg, ccd.CCDData.read(flat_name,unit="adu"), min_value=None, add_keyword=True)
    	print "$$$$$$$$$$  flat corr done $$$$$$$$$$$$$$$", sci_img


    	#the unmasked image is saved
    	ccd.fits_ccddata_writer(fcimg,main_path+sci_folder+'/reduced/'+output_name)
    	#ccd.fits_ccddata_writer(fcimg,main_path+sci_folder+'/reduced/pre_'+output_name)

    	del timg
    	del sdimg
    	del fcimg


    	#print "image %s reduced (trimmed, dark sub, flat corr, masked)" % (day+'/'+sci_folder+'/reduced/'+output_name)
    	print "image %s reduced (trimmed, dark sub, flat corr)" % (day+'/'+sci_folder+'/reduced/'+output_name)

    return(dark_name)



def multi_run_wrapper(args):
#function necessary for the use of pool.map with different arguments
   return red_img(*args)





for j, item in enumerate(folders):

    # the darks are combined daily
    fold=folders[j]

    nn=np.where(list_folder==fold)

    fold_epoch=list_epoch[nn]
    fold_name=list_name[nn]
    fold_chips=list_chips[nn]


    arg_list=[]
#the array with the arguments is generatedm this is necessary to use pool.map
    for i, item in enumerate(fold_name):
        arg_list.append((fold_name[i],fold,fold_chips[i],fold_epoch[i]))



    #the function to reduce the images is run:
    pool = Pool(processes=ncores)

    results = pool.map(multi_run_wrapper,arg_list)

    pool.close()
    pool.join()



elapsed = (time.clock() - start)
print elapsed
