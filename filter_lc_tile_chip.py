import scipy.stats as stat
import numpy as np
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
#from matplotlib import cm
from scipy import interpolate
import scipy.signal as sgl
import pyfits as pf
import os
#import matplotlib.pyplot as plt
from decimal import Decimal
from astropy import stats
from scipy.integrate import quad
import multiprocessing as mp



################################################################################

###################################################################
#parameters for the binning

field='Stripe82' # COSMOS, Stripe82, ElaisS1, XMM_LSS, ECDFS

lc_path='/paula1/QUEST/cand_light_curves_Sep2018/'+field+'/'

agn=np.loadtxt(lc_path+'list.txt',dtype='str').transpose()
agn_list=agn.tolist()

bin_name='clean'

###################################################################
#hacemos el filtro


def filt_curve(agn):

    #for i in range(1,30):
    print agn
    arch=pf.open(lc_path+agn)
    datos=arch[1].data
    head0=arch[0].header
    mag=datos['Q']
    errmag=datos['errQ']
    jd0=datos['JD']
    catalog=datos['catalog']
    chip=datos['chip']
    tile=datos['tile']
    fit_rms=datos['fit_rms']
    chi_red=datos['chi_red']
    num_stars=datos['num_stars']
    m=np.where((mag<30) & (chi_red<=20) & (fit_rms<=0.1) & (num_stars>=20))
    jd=jd0[m]
    mag1=mag[m]
    err1=errmag[m]
    catalog1=catalog[m]
    chip1=chip[m]
    tile1=tile[m]
    fit_rms1=fit_rms[m]
    chi_red1=chi_red[m]
    num_stars1=num_stars[m]

    errmedian1=np.median(err1)
    l=np.where((err1<2.0*errmedian1))
    jd=jd[l]
    mag1=mag1[l]
    err1=err1[l]
    catalog1=catalog1[l]
    chip1=chip1[l]
    tile1=tile1[l]
    fit_rms1=fit_rms1[l]
    chi_red1=chi_red1[l]
    num_stars1=num_stars1[l]

    #se ajusta un polinomio de orden 5.
    if len(jd)>0:
        coefficients = np.polyfit(jd, mag1, 5)
        polynomial = np.poly1d(coefficients)
        pol=polynomial(jd)
        #se calcula el rms
        rm=(mag1-pol)**2
        rms=np.sum(rm)
        rms=rms/(len(jd))
        #se hace el filtro por 3 sigmas
        dist=np.abs(mag1-pol)
        sigma=np.sqrt(rms+err1**2)
        n=np.where((dist/sigma)<=1.5)
        dia=jd[n]
        magf1=mag1[n]
        errf1=err1[n]

        catalog1=catalog1[n]
        chip1=chip1[n]
        tile1=tile1[n]
        fit_rms1=fit_rms1[n]
        chi_red1=chi_red1[n]
        num_stars1=num_stars1[n]

        dif=len(jd0)-len(dia)
        print dif

        fluxf2=10**(-0.4*(magf1+48.6))
        errflux2=0.4*np.log(10)*(10.0**(-0.4*(magf1+48.6)))*errf1


        list_chips=list(set(chip1))
        num_chips=len(list_chips)
        num_tiles=len(set(tile1))


        #save the data for sources with detections in only one chip
        if (num_chips==1):

            #se guardan los fits
            c1=pf.Column(name='JD',format='D',array=dia)
            c2=pf.Column(name='Q',format='D',array=magf1)
            c3=pf.Column(name='errQ',format='D',array=errf1)
            c4=pf.Column(name='fluxQ',format='D',array=fluxf2)
            c5=pf.Column(name='errfluxQ',format='D',array=errflux2)
            c6=pf.Column(name='catalog',format='32A',array=catalog1)
            c7=pf.Column(name='chip',format='3A',array=chip1)
            c8=pf.Column(name='tile',format='E',array=tile1)
            c9=pf.Column(name='fit_rms',format='E',array=fit_rms1)
            c10=pf.Column(name='chi_red',format='E',array=chi_red1)
            c11=pf.Column(name='num_stars',format='E',array=num_stars1)

            coldef=pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11])
            thdu=pf.new_table(coldef)
            thdu.writeto(bin_name+'_onechip_'+agn)
            arch1=pf.open(bin_name+'_onechip_'+agn,mode='update')
            head=arch1[0].header

            head['ALPHA']=head0['ALPHA']
            head['DELTA']=head0['DELTA']
            head['uMAG']=head0['uMAG']
            head['uERR']=   head0['uERR']
            head['gMAG']=   head0['gMAG']
            head['gERR']=   head0['gERR']
            head['rMAG']=   head0['rMAG']
            head['rERR']=   head0['rERR']
            head['iMAG']=   head0['iMAG']
            head['iERR']=   head0['iERR']
            head['zMAG']=   head0['zMAG']
            head['zERR']=   head0['zERR']
            head['DELPOINT']=   dif
            head['FILTER']=     'Q'
            head['T_RANGE']=	dia[-1]-dia[0]
            head['T_length']=	len(dia)
            head['CHIP']= chip1[0]
            head['TILE']= tile1[0]

            arch1.flush()

            cmd='mv '+bin_name+'_onechip_'+agn+'  '+lc_path
            print "num epochs= %d" % (len(dia))
            print "lc length= %f" % (dia[-1]-dia[0])
            print "num chips= %d" % num_chips
            os.system(cmd)

        #save the data for sources with detections in only one chip
        elif (num_chips>1):

            nchip=np.arange(len(chip1))

            if num_chips==2:
                lchip1=len(np.where(chip1==list_chips[0])[0])
                lchip2=len(np.where(chip1==list_chips[1])[0])
                lchiptot=len(chip1)

                fracchip1=float(lchip1)/float(lchiptot)
                fracchip2=float(lchip2)/float(lchiptot)

                if (fracchip1>=0.2) and (fracchip2>=0.2):
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[1]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2<0.2):
                    nchip=np.where((chip1==list_chips[0]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn
                elif (fracchip1<0.2) and (fracchip2>=0.2):
                    nchip=np.where((chip1==list_chips[1]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn

            elif num_chips==3:
                lchip1=len(np.where(chip1==list_chips[0])[0])
                lchip2=len(np.where(chip1==list_chips[1])[0])
                lchip3=len(np.where(chip1==list_chips[2])[0])
                lchiptot=len(chip1)

                fracchip1=float(lchip1)/float(lchiptot)
                fracchip2=float(lchip2)/float(lchiptot)
                fracchip3=float(lchip3)/float(lchiptot)


                if (fracchip1>=0.2) and (fracchip2>=0.2) and (fracchip3>=0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[1]) | (chip1==list_chips[2]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2>=0.2) and (fracchip3<0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[1]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2<0.2) and (fracchip3>=0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[2]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1<0.2) and (fracchip2>=0.2) and (fracchip3>=0.2) :
                    nchip=np.where((chip1==list_chips[1]) | (chip1==list_chips[2]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2<0.2) and (fracchip3<0.2) :
                    nchip=np.where((chip1==list_chips[0]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn
                elif (fracchip1<0.2) and (fracchip2>=0.2) and (fracchip3<0.2) :
                    nchip=np.where((chip1==list_chips[1]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn
                elif (fracchip1<0.2) and (fracchip2<0.2) and (fracchip3>=0.2) :
                    nchip=np.where((chip1==list_chips[2]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn



            elif num_chips==4:
                lchip1=len(np.where(chip1==list_chips[0])[0])
                lchip2=len(np.where(chip1==list_chips[1])[0])
                lchip3=len(np.where(chip1==list_chips[2])[0])
                lchip4=len(np.where(chip1==list_chips[3])[0])
                lchiptot=len(chip1)

                fracchip1=float(lchip1)/float(lchiptot)
                fracchip2=float(lchip2)/float(lchiptot)
                fracchip3=float(lchip3)/float(lchiptot)
                fracchip4=float(lchip4)/float(lchiptot)


                if (fracchip1>=0.2) and (fracchip2>=0.2) and (fracchip3>=0.2) and (fracchip4>=0.2):
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[1]) | (chip1==list_chips[2]) | (chip1==list_chips[3]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2>=0.2) and (fracchip3>=0.2) and (fracchip4<0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[1]) | (chip1==list_chips[2]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2>=0.2) and (fracchip3<0.2) and (fracchip4>=0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[1]) | (chip1==list_chips[3]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2<0.2) and (fracchip3>=0.2) and (fracchip4>=0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[2]) | (chip1==list_chips[3]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1<0.2) and (fracchip2>=0.2) and (fracchip3>=0.2) and (fracchip4>=0.2) :
                    nchip=np.where((chip1==list_chips[1]) | (chip1==list_chips[2]) | (chip1==list_chips[3]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2>=0.2) and (fracchip3<0.2) and (fracchip4<0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[1]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2<0.2) and (fracchip3>=0.2) and (fracchip4<0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[2]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2<0.2) and (fracchip3<0.2) and (fracchip4>=0.2) :
                    nchip=np.where((chip1==list_chips[0]) | (chip1==list_chips[3]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1<0.2) and (fracchip2>=0.2) and (fracchip3>=0.2) and (fracchip4<0.2) :
                    nchip=np.where((chip1==list_chips[1]) | (chip1==list_chips[2]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1<0.2) and (fracchip2>=0.2) and (fracchip3<0.2) and (fracchip4>=0.2) :
                    nchip=np.where((chip1==list_chips[1]) | (chip1==list_chips[3]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1<0.2) and (fracchip2<0.2) and (fracchip3>=0.2) and (fracchip4>=0.2) :
                    nchip=np.where((chip1==list_chips[2]) | (chip1==list_chips[3]))
                    errnew=0.05
                    new_name=bin_name+'_morechip_'+agn
                elif (fracchip1>=0.2) and (fracchip2<0.2) and (fracchip3<0.2) and (fracchip4<0.2) :
                    nchip=np.where((chip1==list_chips[0]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn
                elif (fracchip1<0.2) and (fracchip2>=0.2) and (fracchip3<0.2) and (fracchip4<0.2) :
                    nchip=np.where((chip1==list_chips[1]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn
                elif (fracchip1<0.2) and (fracchip2<0.2) and (fracchip3>=0.2) and (fracchip4<0.2) :
                    nchip=np.where((chip1==list_chips[2]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn
                elif (fracchip1<0.2) and (fracchip2<0.2) and (fracchip3<0.2) and (fracchip4>=0.2) :
                    nchip=np.where((chip1==list_chips[3]))
                    errnew=0.0
                    new_name=bin_name+'_onechip_'+agn

            else:
                nchip=np.arange(0,len(chip1),1)
		errnew=0.05
                new_name=bin_name+'_severalchip_'+agn


            dia=dia[nchip]
            magf1=magf1[nchip]
            errf1=errf1[nchip]

            catalog1=catalog1[nchip]
            chip1=chip1[nchip]
            tile1=tile1[nchip]
            fit_rms1=fit_rms1[nchip]
            chi_red1=chi_red1[nchip]
            num_stars1=num_stars1[nchip]

            fluxf2=10**(-0.4*(magf1+48.6))
            errf1=np.sqrt(errf1**2+(errnew)**2)
            errflux2=0.4*np.log(10)*(10.0**(-0.4*(magf1+48.6)))*errf1

            #se guardan los fits
            c1=pf.Column(name='JD',format='D',array=dia)
            c2=pf.Column(name='Q',format='D',array=magf1)
            c3=pf.Column(name='errQ',format='D',array=errf1)
            c4=pf.Column(name='fluxQ',format='D',array=fluxf2)
            c5=pf.Column(name='errfluxQ',format='D',array=errflux2)
            c6=pf.Column(name='catalog',format='32A',array=catalog1)
            c7=pf.Column(name='chip',format='3A',array=chip1)
            c8=pf.Column(name='tile',format='D',array=tile1)
            c9=pf.Column(name='fit_rms',format='D',array=fit_rms1)
            c10=pf.Column(name='chi_red',format='D',array=chi_red1)
            c11=pf.Column(name='num_stars',format='D',array=num_stars1)

            coldef=pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11])
            thdu=pf.new_table(coldef)
            thdu.writeto(new_name)
            arch1=pf.open(new_name,mode='update')
            head=arch1[0].header

            head['ALPHA']=head0['ALPHA']
            head['DELTA']=head0['DELTA']
            head['uMAG']=head0['uMAG']
            head['uERR']=   head0['uERR']
            head['gMAG']=   head0['gMAG']
            head['gERR']=   head0['gERR']
            head['rMAG']=   head0['rMAG']
            head['rERR']=   head0['rERR']
            head['iMAG']=   head0['iMAG']
            head['iERR']=   head0['iERR']
            head['zMAG']=   head0['zMAG']
            head['zERR']=   head0['zERR']
            head['DELPOINT']=   dif
            head['FILTER']=     'Q'
            head['T_RANGE']=	dia[-1]-dia[0]
            head['T_length']=	len(dia)

            arch1.flush()

            cmd='mv '+new_name+'  '+lc_path
            print "num epochs= %d" % (len(dia))
            print "lc length= %f" % (dia[-1]-dia[0])
            print "num chips= %d" % num_chips
            if len(dia)>10: os.system(cmd)
            else: os.system('rm '+new_name)



    arch.close()




pool = mp.Pool(processes=20)

results = pool.map(filt_curve,agn_list)

pool.close()
pool.join()
