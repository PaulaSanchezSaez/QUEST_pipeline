#code that combine the darks and create a master dark for a period
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

ncores=28 #number of cores used to combine the darks. A high fraction of the total is not recomended, since some process use multiple cores.

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

img_type='d.'

list_darks=sorted(glob.glob(main_path+period+'*/*'+img_type+'*/*'+img_type+'*fits')) #list with the fits images available

list_darks=np.core.defchararray.replace(list_darks,main_path,'')

list_darks=np.array(list_darks)

list_epoch=[x[0:8] for x in list_darks] #list with the epoch for every image

list_epoch=np.array(list_epoch)

list_name=[x[27:51] for x in list_darks] #list with the name for every image

list_name=np.array(list_name)

epochs=sorted(glob.glob(main_path+period+'*')) #array with the available epochs

epochs=np.core.defchararray.replace(epochs,main_path,'')

epochs=np.array(epochs)

print "epochs for which the masterdark will be generated = ", epochs


chips_name=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']

chips_name=np.array(chips_name)



def read_data(chip,list_img_chip):
#function to read the data in the ccdproc format
    data_chip=[]

    for k, item in enumerate(list_img_chip):
        tmdark=ccd.trim_image(ccd.CCDData.read(list_img_chip[k],unit="adu"),fits_section='[5:595,5:2395]')
        data_chip.append(tmdark)
        #data_chip.append(ccd.CCDData.read(list_img_chip[k],unit="adu"))
        del tmdark
    return(data_chip)


def combine_dark(day,chip,exptime,data_chip):
#function to combine the darks for certain chip at certain exposure time
    output_name='Dark_'+exptime+'_'+chip+'.fits'


    if os.path.isfile(main_path+day+'/darks/'+output_name):
        cmd='rm '+main_path+day+'/darks/'+output_name
        os.system(cmd)
        print "file %s deleted" % (main_path+day+'/darks/'+output_name)

    print main_path+day+'/darks/'+output_name
    mdark=ccd.combine(data_chip,output_file=main_path+day+'/darks/'+output_name,method=u'average',mem_limit=1600000000000,sigma_clip=True,sigma_clip_low_thresh=3,sigma_clip_high_thresh=3,sigma_clip_func=np.ma.median,sigma_clip_dev_func=np.ma.std,dtype=np.float64)
    #tmdark=ccd.trim_image(mdark,fits_section='[5:595,5:2395]')
    #ccd.fits_ccddata_writer(mdark,main_path+day+'/darks/'+output_name)
    #cmd='rm '+main_path+day+'/darks/nt_'+output_name
    #os.system(cmd)
    del mdark
    #del tmdark


    #cmd='mv '+output_name+' '+main_path+day+'/darks'
    #os.system(cmd)
    print "masterdark %s created for chip %s for epoch %s" % (exptime,chip,day)

    return(output_name)


def do_comb(chip,list_img_chip,exptime):
#function that runs everithing
    data_chip=read_data(chip,list_img_chip)
    output_name=combine_dark(day,chip,exptime,data_chip)
    del data_chip[:] ; del data_chip
    return(output_name)


def multi_run_wrapper(args):
#function necessary for the use of pool.map with different arguments
   return do_comb(*args)

#the darks are combined considering the day, exptime and chip



for j, item in enumerate(epochs):
    # the darks are combined daily
    day=epochs[j]

    if day in list_epoch:

        if (os.path.exists(main_path+day+'/darks')==False):
            os.system('mkdir '+main_path+day+'/darks')
            print "new dark folder created in %s " % (day)

        nn=np.where(list_epoch==day)

        day_name=list_name[nn]
        day_darks=list_darks[nn]

        #the darks are separated by their exptime

        list_10s=[]
        chip_10s=[]
        list_60s=[]
        chip_60s=[]
        list_180s=[]
        chip_180s=[]

        for i, item in enumerate(day_name):
            a=pf.open(main_path+day_darks[i])
            head=a[0].header
            exptime=head['EXPTIME']
            imgtype=head['IMAGETYP']
            chip_row=head['CHIP-ROW']
            chip_col=head['CHIP-COL']
            if len(str(chip_col+1))==1: chip_col='0'+str(chip_col+1)
            else: chip_col=str(chip_col+1)
            chip=chip_row+chip_col

            if (int(exptime)==10) and (imgtype=='dark'):
                list_10s.append(main_path+day_darks[i])
                chip_10s.append(chip)

            elif (int(exptime)==60) and (imgtype=='dark'):
                list_60s.append(main_path+day_darks[i])
                chip_60s.append(chip)

            elif (int(exptime)==180) and (imgtype=='dark'):
                list_180s.append(main_path+day_darks[i])
                chip_180s.append(chip)

            a.close()
            del a

        list_10s=np.array(list_10s)
        chip_10s=np.array(chip_10s)
        list_60s=np.array(list_60s)
        chip_60s=np.array(chip_60s)
        list_180s=np.array(list_180s)
        chip_180s=np.array(chip_180s)


        #the darks are separated by chip and then they are combined





        #do everithing for 10s

        if len(list_10s)>0:

            pool = Pool(processes=ncores) #56 cores are used

            arg_list=[]
            #the array with the arguments is generatedm this is necessary to use pool.map
            for i, item in enumerate(chips_name):
                list_img_chip=list_10s[np.where(chip_10s==chips_name[i])]
                arg_list.append((chips_name[i],list_img_chip,'10s'))

            results = pool.map(multi_run_wrapper,arg_list)


            arg_list=[]
            del list_10s
            del chip_10s
            del arg_list

            pool.close()
            pool.join()

            print "10s done"

        #do everithing for 60s
        if len(list_60s)>0:

            pool = Pool(processes=ncores) #56 cores are used

            arg_list=[]
            #the array with the arguments is generatedm this is necessary to use pool.map
            for i, item in enumerate(chips_name):
                list_img_chip=list_60s[np.where(chip_60s==chips_name[i])]

                arg_list.append((chips_name[i],list_img_chip,'60s'))

            results = pool.map(multi_run_wrapper,arg_list)

            arg_list=[]
            del list_60s
            del chip_60s
            del arg_list

            pool.close()
            pool.join()

            print "60s done"


        #do everithing for 180s

        if len(list_180s)>0:

            pool = Pool(processes=ncores) #56 cores are used

            arg_list=[]
            #the array with the arguments is generatedm this is necessary to use pool.map
            for i, item in enumerate(chips_name):
                list_img_chip=list_180s[np.where(chip_180s==chips_name[i])]

                arg_list.append((chips_name[i],list_img_chip,'180s'))

            results = pool.map(multi_run_wrapper,arg_list)

            arg_list=[]
            del list_180s
            del chip_180s
            del arg_list

            pool.close()
            pool.join()

            print "180s done"


elapsed = (time.clock() - start)
print elapsed
