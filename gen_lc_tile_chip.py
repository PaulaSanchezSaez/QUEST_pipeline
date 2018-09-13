#generate light curves for objects in a catalog, considering the external flags.
import numpy as np
import pyfits as pf
import math
import os
import time
from scipy.spatial import cKDTree
from astropy import units as u
#from astropy.coordinates import ICRS
from astropy.coordinates import SkyCoord
from astropy.coordinates.matching import match_coordinates_sky as match_coord
import multiprocessing as mp
import pandas as pd

start = time.clock()

field='ECDFS' # COSMOS, Stripe82, ElaisS1, XMM_LSS, ECDFS


cat_path='/paula1/QUEST/catalogs_Sep2018/'+field+'/'
lc_path='/paula1/QUEST/light_curves_Sep2018/'+field+'/'
defa='/paula3/QUEST/scripts/default_files/'
#lc_path='/paula1/QUEST/cand_light_curves_Sep2018/'+field+'/'

dra=0.01 #size of the seccion eliminated from the images in ra
ddec=0.01 #size of the seccion eliminated from the images in dec



###################################################################
#read the photometric catalog (sdss o des)


catq=pf.open(defa+'photometric_catalogs/cat_'+field+'_clean_all_Q.fits')
dat=catq[1].data
raq=dat['ra']
decq=dat['dec']
'''

input_file='/paula1/QUEST/candidates_run_Sep2018/candidates_qso_hp_ws_m16_rand_forest_'+field+'_bin3_May2018.txt'
df = pd.read_csv(input_file)
raq=df['RA']
decq=df['DEC']
'''

if (field=='COSMOS') or (field=='Stripe82') or (field=='XMM_LSS'):
    uu=dat['u']
    uerr=dat['ERR_u']
else:
    uu=np.ones(len(raq))*(-99.0)
    uerr=np.ones(len(raq))*(-99.0)


g=dat['g']
gerr=dat['ERR_g']
r=dat['r']
rerr=dat['ERR_r']
ii=dat['i']
ierr=dat['ERR_i']
z=dat['z']
zerr=dat['ERR_z']

'''

g=df['g']
r=df['r']
ii=df['i']
z=df['z']

gerr=np.ones(len(raq))*(-99.0)
rerr=np.ones(len(raq))*(-99.0)
ierr=np.ones(len(raq))*(-99.0)
zerr=np.ones(len(raq))*(-99.0)
'''



n=int((float(len(raq))/5.0))

'''
raq=raq[0:n]
decq=decq[0:n]
g=g[0:n]
r=r[0:n]
ii=ii[0:n]
uu=uu[0:n]
z=z[0:n]
gerr=gerr[0:n]
rerr=rerr[0:n]
ierr=ierr[0:n]
uerr=uerr[0:n]
zerr=zerr[0:n]
'''

raq=raq[n:2*n]
decq=decq[n:2*n]
g=g[n:2*n]
r=r[n:2*n]
ii=ii[n:2*n]
uu=uu[n:2*n]
z=z[n:2*n]
gerr=gerr[n:2*n]
rerr=rerr[n:2*n]
ierr=ierr[n:2*n]
uerr=uerr[n:2*n]
zerr=zerr[n:2*n]

'''
raq=raq[2*n:3*n]
decq=decq[2*n:3*n]
g=g[2*n:3*n]
r=r[2*n:3*n]
ii=ii[2*n:3*n]
uu=uu[2*n:3*n]
z=z[2*n:3*n]
gerr=gerr[2*n:3*n]
rerr=rerr[2*n:3*n]
ierr=ierr[2*n:3*n]
uerr=uerr[2*n:3*n]
zerr=zerr[2*n:3*n]


raq=raq[3*n:4*n]
decq=decq[3*n:4*n]
g=g[3*n:4*n]
r=r[3*n:4*n]
ii=ii[3*n:4*n]
uu=uu[3*n:4*n]
z=z[3*n:4*n]
gerr=gerr[3*n:4*n]
rerr=rerr[3*n:4*n]
ierr=ierr[3*n:4*n]
uerr=uerr[3*n:4*n]
zerr=zerr[3*n:4*n]




raq=raq[4*n:]
decq=decq[4*n:]
g=g[4*n:]
r=r[4*n:]
ii=ii[4*n:]
uu=uu[4*n:]
z=z[4*n:]
gerr=gerr[4*n:]
rerr=rerr[4*n:]
ierr=ierr[4*n:]
uerr=uerr[4*n:]
zerr=zerr[4*n:]

'''


#convert units
cat_coord = SkyCoord(ra=raq*u.degree, dec=decq*u.degree)

#list for with every objet in the catalog:
#d=np.empty((len(raq),0)).tolist()


#read the list of catalogs



rlist=np.loadtxt(cat_path+'cat_list.txt',dtype='str').transpose()



#namelist=np.core.defchararray.replace(rlist,'.phot.fits','')
#namelist=np.core.defchararray.replace(namelist,cat_path,'')
#nameimg=np.core.defchararray.replace(namelist,'R_','')
#name=np.core.defchararray.rsplit(nameimg,'.')
###################################################################
#iteramos por imagen:

def cat_match(rlist,ncores,nthread):

    ncat=int(len(rlist)/ncores)

    d=np.empty((len(raq),0)).tolist()
    dchips=np.empty((len(raq),0)).tolist()

    if nthread<(ncores-1):
        rlist=rlist[nthread*ncat:(nthread+1)*ncat]
        #name=name[nthread*ncat:(nthread+1)*ncat]

    else:
        rlist=rlist[nthread*ncat:-1]
        #name=name[nthread*ncat:-1]


    largo=len(rlist)
    for j, item in enumerate(rlist):
        if ((".A26." not in rlist[j])  and (".C26." not in rlist[j]) and (".A05." not in rlist[j]) and (".A10." not in rlist[j]) and (".A12." not in rlist[j]) and (".A15." not in rlist[j])  and (".A28." not in rlist[j]) and (".B01." not in rlist[j]) and (".B03." not in rlist[j]) and (".B04." not in rlist[j]) and (".B07." not in rlist[j]) and (".C06." not in rlist[j]) and (".C07." not in rlist[j]) and (".C16." not in rlist[j]) and (".D07." not in rlist[j])):
            arch=pf.open(cat_path+rlist[j])
            datos=arch[1].data
            head=arch[0].header

            alpha=datos['ra']
            delta=datos['dec']
            flags=datos['FLAGS']
            imaflags=datos['IMAFLAGS_ISO']
            nimaflags=datos['NIMAFLAGS_ISO']

            jd=float(head['JD'])
            field=head['Field']
            chip=head['CHIP']
            tile=head['Tile']
            fit_rms=head['fit_rms']
            chi_red=head['chi_red']
            num_stars=head['num_stars']


            if len(alpha)>1:
                min_ra=np.amin(alpha)
                max_ra=np.amax(alpha)
                min_dec=np.amin(delta)
                max_dec=np.amax(delta)


                img_coord = SkyCoord(ra=alpha*u.degree, dec=delta*u.degree)
                ind,ang,dis=match_coord(cat_coord,img_coord,nthneighbor=1)
                ang0=np.array(ang)
                #n=np.where(ang0<0.00008)
                n=np.where(ang0<0.000277778)
                n=n[0] #posicion en el catalogo de los objetos
                if len(n)>0:
                    pos=ind[n] #posicion en la imagen
                    alpha0=alpha[pos]
                    delta0=delta[pos]
                    flags0=flags[pos]
                    imaflags0=imaflags[pos]
                    nimaflags0=nimaflags[pos]
                    mag=datos['Q'][pos]
                    err=datos['errQ'][pos]



                    for i, item in enumerate(pos):
                        if (flags0[i]==0) and (imaflags0[i]==0) and (nimaflags0[i]==0) and (alpha0[i]>(min_ra+dra)) and (alpha0[i]<(max_ra-dra)) and (delta0[i]>(min_dec+ddec)) and (delta0[i]<(max_dec-ddec)):

                            e=[jd,mag[i],err[i]]
                            d[n[i]].append(e)
                            dchips[n[i]].append([rlist[j],chip,tile,fit_rms,chi_red,num_stars])
                            #print d[ind[n][i]]

            print "done image %s    %f     %d" % (rlist[j],jd, (largo-j))
    print len(d)
    return (d,dchips)



#cat_match(rlist,1,0)


pool=mp.Pool(processes=20)
p0=pool.apply_async(cat_match,args=(rlist,20,0))
p1=pool.apply_async(cat_match,args=(rlist,20,1))
p2=pool.apply_async(cat_match,args=(rlist,20,2))
p3=pool.apply_async(cat_match,args=(rlist,20,3))
p4=pool.apply_async(cat_match,args=(rlist,20,4))
p5=pool.apply_async(cat_match,args=(rlist,20,5))
p6=pool.apply_async(cat_match,args=(rlist,20,6))
p7=pool.apply_async(cat_match,args=(rlist,20,7))
p8=pool.apply_async(cat_match,args=(rlist,20,8))
p9=pool.apply_async(cat_match,args=(rlist,20,9))
p10=pool.apply_async(cat_match,args=(rlist,20,10))
p11=pool.apply_async(cat_match,args=(rlist,20,11))
p12=pool.apply_async(cat_match,args=(rlist,20,12))
p13=pool.apply_async(cat_match,args=(rlist,20,13))
p14=pool.apply_async(cat_match,args=(rlist,20,14))
p15=pool.apply_async(cat_match,args=(rlist,20,15))
p16=pool.apply_async(cat_match,args=(rlist,20,16))
p17=pool.apply_async(cat_match,args=(rlist,20,17))
p18=pool.apply_async(cat_match,args=(rlist,20,18))
p19=pool.apply_async(cat_match,args=(rlist,20,19))



'''
pool=mp.Pool(processes=8)
p0=pool.apply_async(cat_match,args=(rlist,8,0))
p1=pool.apply_async(cat_match,args=(rlist,8,1))
p2=pool.apply_async(cat_match,args=(rlist,8,2))
p3=pool.apply_async(cat_match,args=(rlist,8,3))
p4=pool.apply_async(cat_match,args=(rlist,8,4))
p5=pool.apply_async(cat_match,args=(rlist,8,5))
p6=pool.apply_async(cat_match,args=(rlist,8,6))
p7=pool.apply_async(cat_match,args=(rlist,8,7))

'''


d0=p0.get()
d1=p1.get()
d2=p2.get()
d3=p3.get()
d4=p4.get()
d5=p5.get()
d6=p6.get()
d7=p7.get()
d8=p8.get()
d9=p9.get()
d10=p10.get()
d11=p11.get()
d12=p12.get()
d13=p13.get()
d14=p14.get()
d15=p15.get()
d16=p16.get()
d17=p17.get()
d18=p18.get()
d19=p19.get()



d=np.empty((len(raq),0)).tolist()
dchips=np.empty((len(raq),0)).tolist()
nd=[]

for i, item in enumerate(raq):
    d[i]=d0[0][i]+d1[0][i]+d2[0][i]+d3[0][i]+d4[0][i]+d5[0][i]+d6[0][i]+d7[0][i]+d8[0][i]+d9[0][i]+d10[0][i]+d11[0][i]+d12[0][i]+d13[0][i]+d14[0][i]+d15[0][i]+d16[0][i]+d17[0][i]+d18[0][i]+d19[0][i]
    dchips[i]=d0[1][i]+d1[1][i]+d2[1][i]+d3[1][i]+d4[1][i]+d5[1][i]+d6[1][i]+d7[1][i]+d8[1][i]+d9[1][i]+d10[1][i]+d11[1][i]+d12[1][i]+d13[1][i]+d14[1][i]+d15[1][i]+d16[1][i]+d17[1][i]+d18[1][i]+d19[1][i]
    nd.append(len(d[i]))



###################################################################
#construction of the light curves

print "match done"



print "start light curve generation"


def gen_lc(d,dchips,ra,dec,umag,uerr,gmag,gerr,rmag,rerr,imag,ierr,zmag,zerr):
    #function to create every light curve
    f=np.array(d)
    fchips=np.array(dchips)
    #print "arreglo",f

    if f.any():
        f=f[f[:,0].argsort()]
        fchips=fchips[f[:,0].argsort()]
        jdo=f[:,0]
        #print jdo
        if (len(jdo)>2):
            c1=pf.Column(name='JD',format='D',array=jdo)
            c2=pf.Column(name='Q',format='D',array=f[:,1])
            c3=pf.Column(name='errQ',format='D',array=f[:,2])
            c4=pf.Column(name='catalog',format='32A',array=fchips[:,0])
            c5=pf.Column(name='chip',format='3A',array=fchips[:,1])
            c6=pf.Column(name='tile',format='D',array=fchips[:,2])
            c7=pf.Column(name='fit_rms',format='D',array=fchips[:,3])
            c8=pf.Column(name='chi_red',format='D',array=fchips[:,4])
            c9=pf.Column(name='num_stars',format='D',array=fchips[:,5])


            coldef=pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9])
            thdu=pf.new_table(coldef)
            ras="%.6f" % ra
            decs="%.6f" % dec
            thdu.writeto(ras+"_"+decs+"_"+field+".fits")
            arch1=pf.open(ras+"_"+decs+"_"+field+".fits",mode='update')
            head=arch1[0].header
            head['ALPHA']=    ra
            head['DELTA']=    dec

            head['uMAG']=  umag
            head['uERR']=  uerr
            head['gMAG']=  gmag
            head['gERR']=  gerr
            #print "gmag", gmag[j]
            head['rMAG']=   rmag
            head['rERR']=  rerr
            #print "rmag", rmag[j]
            head['iMAG']=   imag
            head['iERR']=  ierr
            #print "imag", imag[j]
            head['zMAG']=   zmag
            head['zERR']=  zerr



            #print "hr", hr[j]
            arch1.flush()
            arch1.close()
            del arch1
            cmd="mv "+ras+"_"+decs+"_"+field+".fits "+lc_path
            os.system(cmd)
            print "done obj = %s \t" % (ras+"_"+decs+"_"+field+".fits")
    return (ra)


def run_gen_lc(args):
#function necessary to use of pool.map when running gen_lc with different arguments
   return gen_lc(*args)




args_run_gen_lc=[]
#the array with the arguments is generatedm this is necessary to use pool.map
for i, item in enumerate(raq):
    args_run_gen_lc.append((d[i],dchips[i],raq[i],decq[i],uu[i],uerr[i],g[i],gerr[i],r[i],rerr[i],ii[i],ierr[i],z[i],zerr[i]))


print len(args_run_gen_lc)
#gen_lc is run
pool = mp.Pool(processes=20)

results = pool.map(run_gen_lc,args_run_gen_lc)

pool.close()
pool.join()









###################################################################

elapsed = (time.clock() - start)
print elapsed
