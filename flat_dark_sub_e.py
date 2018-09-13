#code to substract the darks from the eve flat fields
import time
import numpy as np
import pyfits as pf
import os
import glob
from astropy import units as u
import ccdproc as ccd
from multiprocessing import Pool
import sys, getopt

###############################################################################
#modify these and only these variables:

period='201112' #period for which the list will be generated

ncores=28 #number of cores used to subtract the darks. A high fraction of the total is not recomended, since some process use multiple cores.

main_path='/paula3/QUEST/data/'  #where are the daily image folders located

###############################################################################
# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"p:n:d:")
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



period=str(period) #period for which the list will be generated

ncores=int(ncores) #number of cores used in the reduction process

main_path=str(main_path)

###############################################################################


start = time.clock()

img_type='e.'

list_flats=sorted(glob.glob(main_path+period+'*/*'+img_type+'*/*'+img_type+'*fits')) #list with the fits images available

list_flats=np.core.defchararray.replace(list_flats,main_path,'')

list_flats=np.array(list_flats)

list_epoch=[x[0:8] for x in list_flats] #list with the epoch for every image

list_epoch=np.array(list_epoch)

list_name=[x[27:51] for x in list_flats] #list with the name for every image

list_name=np.array(list_name)

epochs=sorted(glob.glob(main_path+period+'*')) #array with the available epochs

epochs=np.core.defchararray.replace(epochs,main_path,'')

epochs=np.array(epochs)

print "epochs for which the masterdark will be substracted from the flat = ", epochs


chips_name=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']

chips_name=np.array(chips_name)




def subs_dark(flat_dir,flat_name,chip,day):
#function to substract de dark from the flat

    dark_name='Dark_10s_'+chip+'.fits'

    if os.path.isfile(main_path+day+'/darks/'+dark_name):

        if os.path.isfile(main_path+day+'/flats/DS_'+flat_name):
            #print "esta la imagen $$$$$$$"
            os.system('rm  '+main_path+day+'/flats/DS_'+flat_name)
            print "image deleted %s" % (main_path+day+'/flats/DS_'+flat_name)
        #the flat is trimmed
        tflat=ccd.trim_image(ccd.CCDData.read(flat_dir,unit="adu"),fits_section='[5:595,5:2395]')
        #print "$$$$$$$$$$  trim done $$$$$$$$$$$$$$$"
        #the darks is substracted
        dsflat=ccd.subtract_dark(tflat,ccd.CCDData.read(main_path+day+'/darks/'+dark_name,unit="adu"),exposure_time='EXPTIME',exposure_unit=u.second,scale=False)
        #print "$$$$$$$$$$  subs dark done $$$$$$$$$$$$$$$"
        ccd.fits_ccddata_writer(dsflat,main_path+day+'/flats/DS_'+flat_name)

        del tflat
        del dsflat

        print "masterdark %s substracted from flat %s for epoch %s" % (dark_name,flat_name,day)

    return(dark_name)

def multi_run_wrapper(args):
#function necessary for the use of pool.map with different arguments
   return subs_dark(*args)

#the darks are substracted from every flat

for j, item in enumerate(epochs):
    # the darks are combined daily
    day=epochs[j]

    if day in list_epoch:

        if (os.path.exists(main_path+day+'/flats')==False):
            os.system('mkdir '+main_path+day+'/flats')
            print "new flats folder created in %s " % (day)

        #the darks are loaded
        #day_darks=sorted(glob.glob(main_path+day+'/darks/*10s*.fits'))
        #day_darks=np.core.defchararray.replace(list_darks,main_path+day+'/darks/,'')


        nn=np.where(list_epoch==day)

        day_name=list_name[nn]
        day_flats=list_flats[nn]


        list_10s=[]
        chip_10s=[]
        name_10s=[]

        for i, item in enumerate(day_name):
            a=pf.open(main_path+day_flats[i])
            head=a[0].header
            exptime=head['EXPTIME']
            imgtype=head['IMAGETYP']
            chip_row=head['CHIP-ROW']
            chip_col=head['CHIP-COL']
            if len(str(chip_col+1))==1: chip_col='0'+str(chip_col+1)
            else: chip_col=str(chip_col+1)
            chip=chip_row+chip_col

            if (int(exptime)==10) and (imgtype=='pmskyflat'):
                list_10s.append(main_path+day_flats[i])
                chip_10s.append(chip)
                name_10s.append(day_name[i])

            else: print "flat with another exptime"

            a.close()
            del a

        name_10s=np.array(name_10s)
        chip_10s=np.array(chip_10s)
        list_10s=np.array(list_10s)


        pool = Pool(ncores) #56 cores are used


        #do everithing for 10s
        arg_list=[]
        #the array with the arguments is generatedm this is necessary to use pool.map
        for i, item in enumerate(name_10s):
            arg_list.append((list_10s[i],name_10s[i],chip_10s[i],day))
            if (chip_10s[i] in name_10s[i])==False: print "$$$$$ALERTTTTT CHIP IS WRONG$$$$"

        results = pool.map(multi_run_wrapper,arg_list)

        arg_list=[]
        del list_10s
        del name_10s
        del chip_10s
        del arg_list

        pool.close()
        pool.join()

elapsed = (time.clock() - start)
print elapsed
