import scipy.stats as stat
import numpy as np
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
import scipy.signal as sgl
import pyfits as pf
import os
#import matplotlib.pyplot as plt
import multiprocessing as mp



################################################################################

###################################################################
#parameters for the binning

field='Stripe82' # COSMOS, Stripe82, ElaisS1, XMM_LSS, ECDFS

lc_path='/paula1/QUEST/cand_light_curves_Sep2018/'+field+'/'

agn=np.loadtxt(lc_path+'clean_list.txt',dtype='str').transpose()
agn=np.core.defchararray.replace(agn,'clean_','')

agn_list=agn.tolist()

delta=3

bin_name='bin'+str(delta)

###################################################################
#hacemos el filtro

def bin_curve(agn):
    arch=pf.open(lc_path+'clean_'+agn)
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
    m=np.where((mag<30) & (chi_red<=20)  & (fit_rms<=0.1) & (num_stars>=20))
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
    if len(jd)>4:
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

        num_chips=len(set(chip1))
        num_tiles=len(set(tile1))

        #save the data for sources with detections in only one chip
        #if (num_chips==1) & (num_tiles==1):

        minjd=int(np.amin(dia))
        maxjd=int(np.amax(dia)+1)


        jdin=minjd


        jdfinal=[]
        fluxfinal=[]
        errfinal=[]

        while jdin < maxjd:
            dates=dia[np.where((dia>jdin) & (dia<(jdin+delta)))]
            fluxes=fluxf2[np.where((dia>jdin) & (dia<(jdin+delta)))]
            errors=errflux2[np.where((dia>jdin) & (dia<(jdin+delta)))]

            if len(dates)==1:
                jdfinal.append(dates[0])
                fluxfinal.append(fluxes[0])
                errfinal.append(errors[0])

            elif len(dates)>1:
                jdf=np.mean(dates)

                fluxes2=fluxes*1e30
                errors2=errors*1e30

                fluxfb=np.mean(fluxes2)

                err_std=np.sqrt((1.0/(len(fluxes2)*(len(fluxes2)-1)))*np.sum((fluxes2-fluxfb)**2.0,dtype=np.float64),dtype=np.float64)

                dgauss=(1.0/len(fluxes2))*np.sqrt(np.sum(errors2**2,dtype=np.float64),dtype=np.float64)

                errf2=np.amax([err_std,dgauss])

                jdfinal.append(jdf)

                fluxf=fluxfb*1e-30
                errf=errf2*1e-30

                #print err_std,dgauss,errf

                fluxfinal.append(fluxf)
                errfinal.append(errf)

            jdin+=delta

        jdfinal=np.array(jdfinal)
        fluxfinal=np.array(fluxfinal)
        errfinal=np.array(errfinal)
        print "jd fin", len(jdfinal)

        magfinal=-2.5*np.log10(fluxfinal)-48.6
        magerrfinal=(2.5/np.log(10))*(errfinal/fluxfinal)

        #se guardan los fits
        c1=pf.Column(name='JD',format='D',array=jdfinal)
        c2=pf.Column(name='Q',format='D',array=magfinal)
        c3=pf.Column(name='errQ',format='D',array=magerrfinal)
        c4=pf.Column(name='fluxQ',format='D',array=fluxfinal)
        c5=pf.Column(name='errfluxQ',format='D',array=errfinal)


        coldef=pf.ColDefs([c1,c2,c3,c4,c5])
        thdu=pf.new_table(coldef)
        thdu.writeto(bin_name+'_'+agn)
        arch1=pf.open(bin_name+'_'+agn,mode='update')
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
        head['T_RANGE']=	jdfinal[-1]-jdfinal[0]
        head['T_length']=	len(jdfinal)
        head['CHIP']= chip1[0]

        arch1.flush()

        cmd='mv '+bin_name+'_'+agn+'  '+lc_path
        print "num epochs= %d" % (len(jdfinal))
        print "lc length= %f" % (jdfinal[-1]-jdfinal[0])
        if len(jdfinal)>10: os.system(cmd)
        else: os.system('rm '+bin_name+'_'+agn)
        print "%s done " % (agn)

        '''
        elif (num_chips>1) |  (num_tiles>1):

            minjd=int(np.amin(dia))
            maxjd=int(np.amax(dia)+1)


            jdin=minjd


            jdfinal=[]
            fluxfinal=[]
            errfinal=[]

            while jdin < maxjd:
                dates=dia[np.where((dia>jdin) & (dia<(jdin+delta)))]
                fluxes=fluxf2[np.where((dia>jdin) & (dia<(jdin+delta)))]
                errors=errflux2[np.where((dia>jdin) & (dia<(jdin+delta)))]

                if len(dates)==1:
                    jdfinal.append(dates[0])
                    fluxfinal.append(fluxes[0])
                    errfinal.append(errors[0])

                elif len(dates)>1:
                    jdf=np.mean(dates)

                    fluxes2=fluxes*1e30
                    errors2=errors*1e30

                    fluxfb=np.mean(fluxes2)

                    err_std=np.sqrt((1.0/(len(fluxes2)*(len(fluxes2)-1)))*np.sum((fluxes2-fluxfb)**2.0,dtype=np.float64),dtype=np.float64)

                    dgauss=(1.0/len(fluxes2))*np.sqrt(np.sum(errors2**2,dtype=np.float64),dtype=np.float64)

                    errf2=np.amax([err_std,dgauss])

                    jdfinal.append(jdf)

                    fluxf=fluxfb*1e-30
                    errf=errf2*1e-30

                    #print err_std,dgauss,errf

                    fluxfinal.append(fluxf)
                    errfinal.append(errf)

                jdin+=delta

            jdfinal=np.array(jdfinal)
            fluxfinal=np.array(fluxfinal)
            errfinal=np.array(errfinal)
            print "jd fin", len(jdfinal)

            magfinal=-2.5*np.log10(fluxfinal)-48.6
            magerrfinal=(2.5/np.log(10))*(errfinal/fluxfinal)

            magerrfinal=np.sqrt(magerrfinal**2+(0.05)**2)
            errfinal=0.4*np.log(10)*(10.0**(-0.4*(magfinal+48.6)))*magerrfinal

            #se guardan los fits
            c1=pf.Column(name='JD',format='D',array=jdfinal)
            c2=pf.Column(name='Q',format='D',array=magfinal)
            c3=pf.Column(name='errQ',format='D',array=magerrfinal)
            c4=pf.Column(name='fluxQ',format='D',array=fluxfinal)
            c5=pf.Column(name='errfluxQ',format='D',array=errfinal)


            coldef=pf.ColDefs([c1,c2,c3,c4,c5])
            thdu=pf.new_table(coldef)
            thdu.writeto(bin_name+'_morechip_'+agn)
            arch1=pf.open(bin_name+'_morechip_'+agn,mode='update')
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
            head['T_RANGE']=	jdfinal[-1]-jdfinal[0]
            head['T_length']=	len(jdfinal)
            head['CHIP']= chip1[0]
            head['TILE']= tile1[0]

            arch1.flush()

            cmd='mv '+bin_name+'_morechip_'+agn+'  '+lc_path
            print "num epochs= %d" % (len(jdfinal))
            print "lc length= %f" % (jdfinal[-1]-jdfinal[0])
            if len(jdfinal)>10: os.system(cmd)
            else: os.system('rm '+bin_name+'_morechip_'+agn)
            print "%s done " % (agn)

        '''

    arch.close()


pool = mp.Pool(processes=20)

results = pool.map(bin_curve,agn_list)

pool.close()
pool.join()
