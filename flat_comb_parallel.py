#code to generate the masterflats for the period specified
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

period='2011122nd' #period of the month considered

start_epoch=20111216 #initial epoch considered in the period to combine

final_epoch=20111231 #final epoch considered in the period to combine


ncores_rf=15 #number of cores used to read the flat info. 15 or less is recomended.

ncores_cf=14 #number of cores used to combine the flats. It is not recomended to use a high fraction of the cores available. In this case, we use 14 of 64

main_path='/paula3/QUEST/data/'  #where are the daily image folders located

###############################################################################
# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"p:s:f:r:c:d:")
###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-p':
        period=a
    elif o == '-s':
        start_epoch=a
    elif o == '-f':
        final_epoch=a
    elif o == '-r':
        ncores_rf=a
    elif o == '-c':
        ncores_cf=a
    elif o == '-d':
        main_path=a


period=str(period) #period for which the list will be generated

start_epoch=int(start_epoch)  #initial epoch considered in the period to combine

final_epoch=int(final_epoch)  #final epoch considered in the period to combine

ncores_rf=int(ncores_rf) #number of cores used in the reduction process

ncores_cf=int(ncores_cf) #number of cores used in the reduction process

main_path=str(main_path)

###############################################################################


start = time.clock()

print "epochs for which the masterflat will be created = %s  " % (period)

chips_name=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']

#chips_name=['C13','C14']

chips_name=np.array(chips_name)

epoch_list=np.linspace(start_epoch,final_epoch, (final_epoch-start_epoch)+1,dtype=np.int) #list with the epochs in the period


#a new folder where the masterflats will be stored is created
if (os.path.exists(main_path+'MasterFlats')==False):
    os.system('mkdir '+main_path+'MasterFlats')
    print "new MasterFlats folder created in %s " % (main_path)

if (os.path.exists(main_path+'MasterFlats/'+period)==False):
    os.system('mkdir '+main_path+'MasterFlats/'+period)
    print "new MasterFlats/%s folder created in %s " % (period,main_path)


if (os.path.exists(main_path+'Masks')==False):
    os.system('mkdir '+main_path+'Masks')
    print "new Masks folder created in %s " % (main_path)



if (os.path.exists(main_path+'Masks/L')==False):
    os.system('mkdir '+main_path+'Masks/L')
    print "new low counts flats folder created in %s " % (main_path+'Masks')


if (os.path.exists(main_path+'Masks/H')==False):
    os.system('mkdir '+main_path+'Masks/H')
    print "new high counts flats folder created in %s " % (main_path+'Masks')

#the statistics for every flat in the period is calculated

#list1=open(main_path+'MasterFlats/'+period+'/flat_info_'+period+'.txt','w')
#list1.writelines("%s \t %s\t %s \t %s \t %s \n " % ('#period','flats','mean','std','chip'))

#flat_list=[]
#flat_mean=[]
#flat_std=[]
#flat_chip=[]


#for j, item in enumerate(epoch_list):
def make_stats(epoch_numb):
    flat_info=[]
    #day=str(epoch_list[j])
    day=str(epoch_numb)
    print "loading data from epoch %s" % (day)

    if os.path.exists(main_path+day+'/flats'):
        day_flats=sorted(glob.glob(main_path+day+'/flats/DS*fits'))
	if len(day_flats)>0:
            #print day_flats
            day_flats=np.core.defchararray.replace(day_flats,main_path,'')
	    #print day_flats
            for i, item in enumerate(day_flats):
                a=pf.open(main_path+day_flats[i])
                head=a[0].header
                exptime=head['EXPTIME']
                imgtype=head['IMAGETYP']
                chip_row=head['CHIP-ROW']
                chip_col=head['CHIP-COL']

                if len(str(chip_col+1))==1: chip_col='0'+str(chip_col+1)
                else: chip_col=str(chip_col+1)
                chip=chip_row+chip_col

                dat=a[0].data

                dat_mean=np.mean(dat)
                dat_std=np.std(dat)

                #print day,day_flats[i],chip

                flat_info.append((day_flats[i],dat_mean,dat_std,chip))
                #list1.writelines("%s \t %s\t %f \t %f \t %s  \n" % (period,day_flats[i],dat_mean,dat_std,chip))
                #flat_list.append(day_flats[i])
                #flat_mean.append(dat_mean)
                #flat_std.append(dat_std)
                #flat_chip.append(chip)

		if ((500<dat_mean) and (dat_mean<600)):
                    cm='cp '+main_path+day_flats[i]+'   '+main_path+'Masks/L'
		    os.system(cm)

		elif ((2500<dat_mean) and (dat_mean<3500)):
                    cm='cp '+main_path+day_flats[i]+'   '+main_path+'Masks/H'
		    os.system(cm)

                a.close()
                del a

    print "end data from epoch %s" % (day)
    return(flat_info)


#pool = Pool(processes=((final_epoch-start_epoch)+1))
pool = Pool(processes=ncores_rf)

results = pool.map(make_stats,epoch_list.tolist())

pool.close()
pool.join()

flat_info=[]

for i in xrange((final_epoch-start_epoch)+1):

    flat_info+=results[i]

flat_info=np.array(flat_info)

print flat_info
print flat_info.shape



flat_list=flat_info[:,0]
flat_mean=flat_info[:,1].astype(np.float)
flat_std=flat_info[:,2].astype(np.float)
flat_chip=flat_info[:,3]

np.savez(main_path+'MasterFlats/'+period+'/flat_info_'+period+'.npz', flat_list, flat_mean,flat_std,flat_chip)

print flat_mean
print flat_chip

#we combine the flats with mean counts above the avarage

print "start masterflat generation"

def read_data(chip,list_img_chip):
#function to read the data in the ccdproc format
    data_chip=[]

    print "reading chip=%s" % (chip)

    for k, item in enumerate(list_img_chip):
        img=ccd.CCDData.read(main_path+list_img_chip[k],unit="adu")
        data_chip.append(img)
        #data_chip.append(ccd.CCDData.read(list_img_chip[k],unit="adu"))
        del img
    return(data_chip)


def flat_combine(chip,chip_flats,chip_means):
    #function to combine the flat for a certain chip

    output_name='Flat_'+period+'_'+chip+'.fits'

    #nn=np.where(flat_chip==chip)
    #chip_flats=flat_list[nn]
    #chip_means=flat_mean[nn]

    if os.path.isfile(main_path+'MasterFlats/'+period+'/'+output_name):
        os.system('rm  '+main_path+'MasterFlats/'+period+'/'+output_name)
        print "image deleted %s" % (main_path+'MasterFlats/'+period+'/'+output_name)

    mean_counts=np.mean(chip_means)

    mm=np.where((chip_means>mean_counts) & (chip_means<15000))
    flat_selected=chip_flats[mm]

    print "###########",output_name,len(flat_selected), mean_counts

    if len(flat_selected)>0:

        data_chip=read_data(chip,flat_selected)

        print "combining chip=%s" % (chip)

        mflat=ccd.combine(data_chip,output_file=main_path+'MasterFlats/'+period+'/'+output_name,method=u'median',mem_limit=1600000000000,sigma_clip=True,sigma_clip_low_thresh=3,sigma_clip_high_thresh=3,sigma_clip_func=np.ma.median,sigma_clip_dev_func=np.ma.std,dtype=np.float64)

        del mflat

        print "masterflat %s created for chip %s for epoch %s" % (output_name,chip,period)

    else:
        print "chip = %s without normal masterdark " % (chip)

        a=pf.open(main_path+chip_flats[0])
        head=a[0].header
        dat=a[0].data
        fake_data=np.ones(dat.shape)
        pf.writeto(main_path+'MasterFlats/'+period+'/'+output_name,data=fake_data,header=head)

        print "fake masterflat %s created for chip %s for epoch %s" % (output_name,chip,period)

    return(output_name)


def multi_run_wrapper(args):
#function necessary for the use of pool.map with different arguments
   return flat_combine(*args)

#the masterflat are generated


arg_list=[]
#the array with the arguments is generatedm this is necessary to use pool.map
for i, item in enumerate(chips_name):
    list_img_chip=flat_list[np.where(flat_chip==chips_name[i])]
    list_mean_chip=flat_mean[np.where(flat_chip==chips_name[i])]
    arg_list.append((chips_name[i],list_img_chip,list_mean_chip))




pool = Pool(processes=ncores_cf)

results = pool.map(multi_run_wrapper,arg_list)

pool.close()
pool.join()

elapsed = (time.clock() - start)
print elapsed
