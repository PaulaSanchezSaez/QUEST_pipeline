#code to generate the masks for the period specified
import time
import numpy as np
import pyfits as pf
import os
import glob
from astropy import units as u
import ccdproc as ccd
from multiprocessing import Pool
from scipy import ndimage
import sys, getopt



###############################################################################
#modify these and only these variables:

period='2015' #period considered

start_epoch=20150101 #initial epoch considered in the period to combine

final_epoch=20151231 #final epoch considered in the period to combine


ncores=14 #number of cores used to combine the flats and to generate the masks. It is not recomended to use a high fraction of the cores available. In this case, we use 14 of 64

main_path='/paula3/QUEST/data/'  #where are the daily image folders located

###############################################################################
# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"p:s:f:n:d:")
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
    elif o == '-n':
        ncores=a
    elif o == '-d':
        main_path=a




period=str(period) #period considered

start_epoch=int(start_epoch) #initial epoch considered in the period to combine

final_epoch=int(final_epoch) #final epoch considered in the period to combine


ncores=int(ncores) #number of cores used to combine the flats and to generate the masks. It is not recomended to use a high fraction of the cores available. In this case, we use 14 of 64

main_path=str(main_path)
###############################################################################

start = time.time()

print "epoch for which the Masks will be created = %s  " % (period)



#chips_name=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']



chips_name=['B01','B07','B22','B24','B24','B26','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']


#chips_name=np.array(chips_name)

epoch_list=np.linspace(start_epoch,final_epoch, (final_epoch-start_epoch)+1,dtype=np.int) #list with the epochs in the period


#a new folder where the masterflats will be stored is created
if (os.path.exists(main_path+'Masks')==False):
    os.system('mkdir '+main_path+'Masks')
    print "new Masks folder created in %s " % (main_path)


if (os.path.exists(main_path+'Masks/'+period)==False):
    os.system('mkdir '+main_path+'Masks/'+period)
    print "new Masks/%s folder created in %s " % (period,main_path)




#we combine the darks with mean counts above the avarage

print "start Masks generation"

def read_data(chip,list_img_chip):
#function to read the data in the ccdproc format
    data_chip=[]
    print "reading chip=%s" % (chip)

    if len(list_img_chip)>30:
        for k, item in enumerate(list_img_chip):
            img=ccd.CCDData.read(list_img_chip[k],unit="adu")
            img_data=img.data
            img_mean=np.mean(img_data)
            if (img_mean > 3000.) or (img_mean < 550.):
                data_chip.append(img)
                #data_chip.append(ccd.CCDData.read(list_img_chip[k],unit="adu"))
            del img


    else:
        for k, item in enumerate(list_img_chip):
            img=ccd.CCDData.read(list_img_chip[k],unit="adu")
            data_chip.append(img)
            #data_chip.append(ccd.CCDData.read(list_img_chip[k],unit="adu"))
            del img
    return(data_chip)


def gen_mask(chip):
    #function to combine the flat for a certain chip, dividing by the counts number

    output_name_L='Flat_'+period+'_Low_'+chip+'.fits'
    output_name_H='Flat_'+period+'_High_'+chip+'.fits'
    output_name='Mask_'+period+'_'+chip+'.fits'
    output_name_vmap='VarMap_'+period+'_'+chip+'.fits'

    if os.path.isfile(main_path+'Masks/'+period+'/'+output_name_L):
        #print "esta la imagen $$$$$$$"
        os.system('rm  '+main_path+'Masks/'+period+'/'+output_name_L)
        print "image deleted %s" % (main_path+'Masks/'+period+'/'+output_name_L)


    if os.path.isfile(main_path+'Masks/'+period+'/'+output_name_H):
        #print "esta la imagen $$$$$$$"
        os.system('rm  '+main_path+'Masks/'+period+'/'+output_name_H)
        print "image deleted %s" % (main_path+'Masks/'+period+'/'+output_name_H)


    if os.path.isfile(main_path+'Masks/'+period+'/'+output_name):
        #print "esta la imagen $$$$$$$"
        os.system('rm  '+main_path+'Masks/'+period+'/'+output_name)
        print "image deleted %s" % (main_path+'Masks/'+period+'/'+output_name)


    if os.path.isfile(main_path+'Masks/'+period+'/'+output_name_vmap):
        #print "esta la imagen $$$$$$$"
        os.system('rm  '+main_path+'Masks/'+period+'/'+output_name_vmap)
        print "image deleted %s" % (main_path+'Masks/'+period+'/'+output_name_vmap)



    low_list=sorted(glob.glob(main_path+'Masks/L/*'+chip+'*fits'))
    
    if len(low_list)>0:
        low_list=np.array(low_list)
        l_list=np.core.defchararray.replace(low_list,main_path+'Masks/L/DS_','')

        low_epochs=[x[0:8] for x in l_list] #list with the epoch for every image
        low_epochs=np.array(low_epochs).astype(int)

        low_list=low_list[np.where((low_epochs>=start_epoch) & (low_epochs<=final_epoch) )]

    print "combining %d flats with low counts number for chip=%s" % (len(low_list),chip)

    high_list=sorted(glob.glob(main_path+'Masks/H/*'+chip+'*fits'))
    
    if len(high_list)>0:
        high_list=np.array(high_list)
        h_list=np.core.defchararray.replace(high_list,main_path+'Masks/H/DS_','')

        high_epochs=[x[0:8] for x in h_list] #list with the epoch for every image
        high_epochs=np.array(high_epochs).astype(int)

        high_list=high_list[np.where((high_epochs>=start_epoch) & (high_epochs<=final_epoch) )]

    print "combining %d flats with high counts number for chip=%s" % (len(high_list),chip)



    #combining low counts flats
    if len(low_list)>0:

        data_chip=read_data(chip,low_list)

        print "combining low counts flats for chip=%s" % (chip)

        lflat=ccd.combine(data_chip,output_file=main_path+'Masks/'+period+'/'+output_name_L,method=u'median',mem_limit=1600000000000,sigma_clip=True,sigma_clip_low_thresh=3,sigma_clip_high_thresh=3,sigma_clip_func=np.ma.median,sigma_clip_dev_func=np.ma.std,dtype=np.float64)

        del lflat
        del data_chip

        print "low counts flat %s created for chip %s for epoch %s" % (output_name_L,chip,period)



    #combining high counts flats
    if len(high_list)>0:

        data_chip=read_data(chip,high_list)

        print "combining high counts flats for chip=%s" % (chip)

        hflat=ccd.combine(data_chip,output_file=main_path+'Masks/'+period+'/'+output_name_H,method=u'median',mem_limit=1600000000000,sigma_clip=True,sigma_clip_low_thresh=3,sigma_clip_high_thresh=3,sigma_clip_func=np.ma.median,sigma_clip_dev_func=np.ma.std,dtype=np.float64)

        del hflat
        del data_chip

        print "high counts flat %s created for chip %s for epoch %s" % (output_name_H,chip,period)


    if (len(high_list)>0) and (len(low_list)>0):

        h=pf.open(main_path+'Masks/'+period+'/'+output_name_H)
        head_h=h[0].header
        dat_h=h[0].data
        mean_h=np.mean(dat_h)
        norm_h=dat_h/mean_h

        l=pf.open(main_path+'Masks/'+period+'/'+output_name_L)
        head_l=l[0].header
        dat_l=l[0].data
        mean_l=np.mean(dat_l)
        norm_l=dat_l/mean_l


        rflat=norm_l/norm_h

        #pf.writeto(main_path+'Masks/'+period+'/rate_'+output_name,data=rflat,header=head_h)

        #newdata, mask=ccd.cosmicray_median(rflat, thresh=5, mbox=11, gbox=0, rbox=0)
        marr=ndimage.median_filter(rflat,size=(40,200))
        rarr=(rflat-marr)
        m=np.where(np.abs(rarr)>0.05)

        aux=np.zeros(rflat.shape)
        aux[m]=1
        aux=np.int32(aux)


        print aux

        pf.writeto(main_path+'Masks/'+period+'/'+output_name,data=aux,header=head_h)
        arch1=pf.open(main_path+'Masks/'+period+'/'+output_name,mode='update')
        head=arch1[0].header
        head['IMAGETYP']=    'mask'
        #head.update('IMAGETYP',    'mask')
        head['Fake_mask']=    'no'
        arch1.flush()
        del arch1

        print "Mask %s created for chip %s for epoch %s" % (output_name,chip,period)


        aux2=np.ones(rflat.shape)
        aux2[m]=1e30
        uax2=aux2.astype(int)

        pf.writeto(main_path+'Masks/'+period+'/'+output_name_vmap,data=aux2,header=head_h)
        arch1=pf.open(main_path+'Masks/'+period+'/'+output_name_vmap,mode='update')
        head=arch1[0].header
        head['IMAGETYP']=    'varmap'
        #head.update('IMAGETYP',    'mask')
        head['Fake_vmap']=    'no'
        arch1.flush()

        print "VarMap %s created for chip %s for epoch %s" % (output_name_vmap,chip,period)

    '''
    else:

        #aa=pf.open(main_path+chip_flats[0])
        #head=aa[0].header
        #dat=aa[0].data
        fake_data=np.zeros(dat.shape)
        pf.writeto(main_path+'Masks/'+period+'/'+output_name,data=fake_data,header=head)

        arch1=pf.open(main_path+'Masks/'+period+'/'+output_name,mode='update')
        head=arch1[0].header
        head['IMAGETYP']=    'mask'
        #head.update('IMAGETYP',    'mask')
        head['Fake_mask']=    'yes'
        arch1.flush()
        del arch1

        print "$$$$$$$$$$$$$$$ fake mask %s created for chip %s for epoch %s  $$$$$$$$$$$$$$$$$$$" % (output_name,chip,period)



        fake_data=np.ones(dat.shape)
        pf.writeto(main_path+'Masks/'+period+'/'+output_name_vmap,data=fake_data,header=head)

        arch1=pf.open(main_path+'Masks/'+period+'/'+output_name_vmap,mode='update')
        head=arch1[0].header
        head['IMAGETYP']=    'varmap'
        #head.update('IMAGETYP',    'mask')
        head['Fake_vmap']=    'yes'
        arch1.flush()
        del arch1

        print "$$$$$$$$$$$$$$$ fake varmap %s created for chip %s for epoch %s  $$$$$$$$$$$$$$$$$$$" % (output_name_vmap,chip,period)
 
    '''

    if os.path.isfile(main_path+'Masks/'+period+'/'+output_name_L):
        #print "esta la imagen $$$$$$$"
        os.system('rm  '+main_path+'Masks/'+period+'/'+output_name_L)
        print "image deleted %s" % (main_path+'Masks/'+period+'/'+output_name_L)


    if os.path.isfile(main_path+'Masks/'+period+'/'+output_name_H):
        #print "esta la imagen $$$$$$$"
        os.system('rm  '+main_path+'Masks/'+period+'/'+output_name_H)
        print "image deleted %s" % (main_path+'Masks/'+period+'/'+output_name_H)


    return(output_name)



#the masterflat are generated


pool = Pool(processes=ncores)

results = pool.map(gen_mask,chips_name)

pool.close()
pool.join()




#gen_mask('A17')


elapsed = (time.time() - start)
print elapsed
