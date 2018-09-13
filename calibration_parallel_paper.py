import numpy as np
import pyfits as pf
import os
import glob
import time
from multiprocessing import Pool
import sys, getopt
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
from astropy.coordinates.matching import match_coordinates_sky as match_coord


###############################################################################
#modify these and only these variables:

period='201104' #period for which the list will be generated

ncores=28 #number of cores used to combine the darks. A high fraction of the total is not recomended, since some process use multiple cores.

main_path='/paula3/QUEST/data/'  #where are the daily image folders located

defa='/paula3/QUEST/scripts/default_files/'


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
        main_path=a


period=str(period) #period for which the list will be generated

ncores=int(ncores) #number of cores used in the reduction process

main_path=str(main_path)

defa=str(defa)


###############################################################################


start = time.clock()


list_gain=[5.405,5.128,5.988,5.154,5.0,8.474,5.263,6.211,4.524,5.0,5.376,0.0,5.747,5.102,0.0,4.926,5.347,5.714,0.0,4.926,5.347,5.524,5.154,5.586,4.608,15.625,5.524,5.0,0.0,5.291,4.739,5.617,6.134,5.434,0.0,5.555,5.235,5.952,5.988,5.128,4.878,5.780,5.882,5.649,5.555,6.134,5.347,6.756,5.405,4.950,5.347,5.076,6.369,5.917,4.739,5.617,5.154,6.410,6.060,5.617,5.050,5.0,0.0,4.444,5.747,5.128,5.376,5.681,5.102,5.347,5.747,5.0,5.813,5.649,5.917,5.714,8.547,6.993,6.993,5.555,5.988,8.849,4.830,5.747,6.024,6.622,6.289,6.134,4.807,9.345,0.0,7.692,5.617,5.494,5.847,5.988,5.376,5.050,5.988,5.524,4.807,4.950,5.649,5.263,4.672,5.524,5.154,5.181,5.780,5.813,7.751,5.405]
list_gain=np.array(list_gain)

chips_name=['A01','A02','A03','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28','D01','D02','D03','D04','D05','D06','D07','D08','D09','D10','D11','D12','D13','D14','D15','D16','D17','D18','D19','D20','D21','D22','D23','D24','D25','D26','D27','D28']
chips_name=np.array(chips_name)

#info of the available science images

sci=np.loadtxt(defa+'datalogs/sci_images_'+period+'.txt',dtype='str',comments='>>').transpose()
sci_type=sci[2]
sci_jd=sci[6]
sci_img=sci[8]
#tilee=sci[12]
sci_field=sci[12]

tilee=np.ones(len(sci_field),dtype=int)

#sci_ep=[x[0:8] for x in sci_img]
#sci_name=[x[8:15] for x in sci_img]

#we change the name of the fileds to a more useful name


tilee[np.where(sci_field=='COSMOS1')]=1
tilee[np.where(sci_field=='COSMOS2')]=2
tilee[np.where(sci_field=='Stripe82_1')]=1
tilee[np.where(sci_field=='Stripe82_2')]=2
tilee[np.where(sci_field=='Elais_S1_1')]=1
tilee[np.where(sci_field=='Elais_S1_2')]=2
tilee[np.where(sci_field=='XMM_LSS1_1')]=1
tilee[np.where(sci_field=='XMM_LSS1_2')]=2
tilee[np.where(sci_field=='XMM_LSS2_1')]=3
tilee[np.where(sci_field=='XMM_LSS2_2')]=4
tilee[np.where(sci_field=='ECDF_S_1')]=1
tilee[np.where(sci_field=='ECDF_S_2')]=2
tilee[np.where(sci_field=='ECDF1')]=1
tilee[np.where(sci_field=='ECDF2')]=2

sci_field[np.where(sci_field=='COSMOS1')]='COSMOS'
sci_field[np.where(sci_field=='COSMOS2')]='COSMOS'
sci_field[np.where(sci_field=='Stripe82_1')]='Stripe82'
sci_field[np.where(sci_field=='Stripe82_2')]='Stripe82'
sci_field[np.where(sci_field=='Elais_S1_1')]='ElaisS1'
sci_field[np.where(sci_field=='Elais_S1_2')]='ElaisS1'
sci_field[np.where(sci_field=='XMM_LSS1_1')]='XMM_LSS'
sci_field[np.where(sci_field=='XMM_LSS1_2')]='XMM_LSS'
sci_field[np.where(sci_field=='XMM_LSS2_1')]='XMM_LSS'
sci_field[np.where(sci_field=='XMM_LSS2_2')]='XMM_LSS'
sci_field[np.where(sci_field=='ECDF_S_1')]='ECDFS'
sci_field[np.where(sci_field=='ECDF_S_2')]='ECDFS'
sci_field[np.where(sci_field=='ECDF1')]='ECDFS'
sci_field[np.where(sci_field=='ECDF2')]='ECDFS'




img_type='s.'

list_catalogs=sorted(glob.glob(main_path+period+'*/*'+img_type+'*/reduced/*.cat.fits'))

namelist=np.core.defchararray.replace(list_catalogs,'.cat.fits','')

list_catalogs=np.core.defchararray.replace(list_catalogs,main_path,'')

list_catalogs=np.array(list_catalogs)



list_epoch=[x[0:8] for x in list_catalogs] #list with the epoch for every image

list_epoch=np.array(list_epoch)

list_dir=[x[0:34] for x in list_catalogs]

list_dir=np.array(list_dir)

list_name=[x[35:65] for x in list_catalogs] #list with the name for every image

list_name=np.array(list_name)

list_chip=[x[18:21] for x in list_name]

list_chip=np.array(list_chip)

list_img=[x[2:17] for x in list_name] # only contains the name, without the chip and .cat.fits

list_img=np.array(list_img)


############################################################################
#function which run the calibration




def run_calib(image,jd,field,tile,chip):
    print "calibrating image ",image

    if os.path.isfile(image+'.phot.fits'): os.system('rm '+image+'.phot.fits')

    arch=pf.open(image+'.cat.fits',ignore_missing_end=True)
    head=arch[0].header
    if head['SEXNFIN']>0:

        dat=arch[1].data
        alpha=dat['ALPHA_J2000']
        delta=dat['DELTA_J2000']
        flags=dat['FLAGS']
        flag_ima=dat['IMAFLAGS_ISO']
        flag_nima=dat['NIMAFLAGS_ISO']
        flagsf=flags
        clas=dat['CLASS_STAR']
        fwhm=dat['FWHM_IMAGE']
        maginst=dat['MAG_APER'][:,4]
        magQi=maginst
        #print np.amin(maginst),np.amax(maginst)
        errinst=dat['MAGERR_APER'][:,4]
        errmagQi=errinst
        img_coord = SkyCoord(ra=alpha*u.degree, dec=delta*u.degree) 
        #img_coord=ICRS(alpha,delta,unit=(u.degree,u.degree))
        ind,ang,dis=match_coord(cat_coord,img_coord,nthneighbor=1)
        #ang=ang*u.degree
        ang0=np.array(ang)
        n=np.where((ang0<0.000277778))
        n=n[0]
        #print ind
        print "min distance= ",np.amin(ang0)
        if len(n)>0:
            pos=ind[n]
            flags=flags[pos]
            maginst=maginst[pos]
            errinst=errinst[pos]
            q_mag=qmag[n]
            q_err=qerr[n]
            fn=np.where((flags==0) & (q_mag>15.5) & (q_mag<20) & (q_err<0.2) & (maginst<20) & (maginst>15.5) & (errinst<0.2))
            #fn=np.where((flags==0) & (q_mag>15) & (q_mag<18) & (q_err<0.2) & (maginst<19) & (maginst>14) & (errinst<0.3))
            #fn=np.where((flags<=2) & (q_mag<25) & (q_err<0.3) & (maginst<17.5) & (maginst>14) & (errinst<0.3))
            maginst=maginst[fn]
            magtest=maginst

            if (len(maginst)>=15):
                errinst=errinst[fn]
                #print "err_inst",errinst
                q_mag=(qmag[n])[fn]
                q_err=(qerr[n])[fn]
                #print "q_err", q_err
                sigma=np.sqrt(q_err**2+errinst**2)
                errf=sigma**2
                #print "sigma",sigma
                diff=maginst-q_mag

                #print "diff",diff
                coefficients = np.polyfit(maginst, diff, 1)
                polynomial = np.poly1d(coefficients)
                fit=polynomial(maginst)
                #print "fit",fit
                res1=np.abs(diff-fit)
                rms=(np.sum((diff-fit)**2))/(len(maginst)-1)
                sigma=np.sqrt(sigma**2+rms)
                l1=np.where((res1/sigma)<=2)
                r1=np.where((res1/sigma)>2)
                print "first clean",len(l1[0])

                #print "res1",res1
                diff2=diff
                diff=diff[l1]
                diff_rej=diff2[r1]
                sigma=sigma[l1]
                errf=errf[l1]
                maginst=maginst[l1]
                q_mag=q_mag[l1]
                fita=fit[l1]

                '''
                plt.plot(magtest,diff2,'ro')
                plt.plot(maginst,diff,'k*')
                plt.plot(magtest,fit,'b-')
                plt.show()
                '''
                if len(diff)>=15:
                    nit=0
                    while (len(diff)>=15) and (nit<=3):
                        if len(diff_rej)>0:

                            coefficients = np.polyfit(maginst, diff, 1)
                            polynomial = np.poly1d(coefficients)
                            fit=polynomial(maginst)
                            res2=np.abs(diff-fit)
                            rms=(np.sum((diff-fit)**2))/(len(maginst)-1)
                            sigma=np.sqrt(sigma**2+rms)
                            l2=np.where((res2/sigma)<=2)
                            r2=np.where((res2/sigma)>2)
                            diff2=diff
                            diff=diff[l2]
                            diff_rej=diff2[r2]
                            sigma=sigma[l2]
                            errf=errf[l2]
                            maginst=maginst[l2]
                            q_mag=q_mag[l2]
                            nit=nit+1
                        else: break

                    coefficients = np.polyfit(maginst, diff, 1)
                    fit=polynomial(maginst)

                    rms=(np.sum((diff-fit)**2))/(len(maginst)-1)

                    chi_red=(np.sum(((diff-fit)**2)/errf))/(len(maginst)-2)
                    print np.sqrt(rms), chi_red, len(maginst)

                    '''
                    plt.plot(maginst,diff,'ro')
                    plt.plot(maginst,fit,'b-')
                    plt.show()
                    '''

                    final_pol=polynomial
                    final_fit=final_pol(magQi)
                    Q_cal=magQi-final_fit
                    Q_err=np.sqrt(errmagQi**2+rms)

                    #writing the new catalogs:
                    c1=pf.Column(name='ra',format='D',array=alpha)
                    c2=pf.Column(name='dec',format='D',array=delta)
                    c3=pf.Column(name='Q',format='D',array=Q_cal)
                    c4=pf.Column(name='errQ',format='D',array=Q_err)
                    c5=pf.Column(name='FLAGS',format='D',array=flagsf)
                    c6=pf.Column(name='CLASS_STAR',format='D',array=clas)
                    c7=pf.Column(name='FWHM_IMAGE',format='D',array=fwhm)
                    c8=pf.Column(name='Q_inst',format='D',array=magQi)
                    c9=pf.Column(name='IMAFLAGS_ISO',format='D',array=flag_ima)
                    c10=pf.Column(name='NIMAFLAGS_ISO',format='D',array=flag_nima)

                    coldef=pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10])
                    thdu=pf.new_table(coldef)
                    thdu.writeto(image+'.phot.fits')
                    arch1=pf.open(image+'.phot.fits',mode='update')
                    head=arch1[0].header
                    head['JD']=    jd
                    head['Field']=    field
                    head['CHIP']= chip
                    head['Tile']= tile
                    head['fit_rms']= np.sqrt(rms)
                    head['chi_red']= chi_red
                    head['num_stars']= len(maginst)
                    arch1.flush()
                    del arch1
                    print 'catalog generated ###########', image+'.phot.fits', chip
                    cmd='rsync --progress '+image+'.phot.fits  /paula1/QUEST/catalogs_Sep2018/'+field+'/'
		    os.system(cmd)
                    print '######## catalog ', image+'.phot.fits', chip, '  copied to /paula1/QUEST/catalogs_Sep2018/'+field+'/' 

    return(image)



def multi_run_wrapper(args):
#function necessary for the use of pool.map with different arguments
   return run_calib(*args)
############################################################################


#function to run the calibration in every field separately
def calib_field(field):

    #read the list of the catalogs of the images

    nf=np.where((field==sci_field))

    sci_img_field=sci_img[nf[0]]
    sci_jd_field=sci_jd[nf[0]]
    tile_field=tilee[nf[0]]
    #print tile_field
    #print sci_jd_field
    #print sci_img_field

    global cat_coord
    global raq
    global decq
    global qmag
    global qerr


    if len(sci_img_field)>0:



        #read the catalog
        catq=pf.open(defa+'photometric_catalogs/cat_'+field+'_clean_stars_Q.fits')
        dat=catq[1].data
        raq=dat['ra']
        decq=dat['dec']
        qmag=dat['Q']
        qerr=dat['ERR_Q']

        #convert units
        #cat_coord=ICRS(raq,decq,unit=(u.degree,u.degree))
        cat_coord = SkyCoord(ra=raq*u.degree, dec=decq*u.degree)
        #print cat_coord
        #print cat_coord
        #SkyCoord(ICRS, ra=ra, dec=dec




        for i, item in enumerate(sci_img_field):
            #print sci_img_field[i]
            namelist_img= [s for s in namelist if sci_img_field[i] in s]
            #nimg=np.where((sci_img_field[i] in namelist))
            #namelist_img=namelist[nimg] #list with the names of all the chip images for a certain image
            #print namelist_img
            if len(namelist_img)>0:
                #we run the calibration for every chip image
                #namelist_img=namelist_img.tolist()
                '''
                for j, item in enumerate(namelist_img):
                    run_calib(namelist_img[j])

                '''
                jd=sci_jd_field[i]
                tile1=tile_field[i]

                arg_list=[]
                short_namelist_img=np.core.defchararray.replace(namelist_img,main_path,'')

                for j, item in enumerate(namelist_img):
                    chipp=short_namelist_img[j][53:56]
                    #print namelist_img[j]
                    #print chipp
                    arg_list.append((namelist_img[j],jd,field,tile1,chipp))

                #we run the calibration using miltiprocessing
                print ncores
                pool = Pool(ncores) #ncores cores are used

                pool.map(multi_run_wrapper,arg_list)

                pool.close()
                pool.join()





#we run the code for every field
#fields=['COSMOS','Stripe82','ElaisS1','XMM_LSS','ECDFS']
#fields=['Stripe82','XMM_LSS']
#fields=['ElaisS1','ECDFS']
#fields=['ElaisS1','ECDFS','COSMOS','Stripe82','XMM_LSS']

fields=['XMM_LSS','Stripe82','COSMOS']

raq=0
decq=0
qmag=0
qerr=0

#convert units
cat_coord=0


for i,  item in enumerate(fields):
    print "calibrating field=%s" % (fields[i])
    calib_field(fields[i])



elapsed = (time.clock() - start)
print elapsed
