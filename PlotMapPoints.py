import healpy as hp
import sys
import matplotlib.pyplot as plt
import numpy as np

def ga2equ(l,b):
    """
    Convert Galactic to Equatorial coordinates (J2000.0)
    (use at own risk)
    
    Input: [l,b] in decimal degrees
    Returns: [ra,dec] in decimal degrees
    
    Source: 
    - Book: "Practical astronomy with your calculator" (Peter Duffett-Smith)
    - Wikipedia "Galactic coordinates"
    
    Tests (examples given on the Wikipedia page):
    >>> ga2equ([0.0, 0.0]).round(3)
    array([ 266.405,  -28.936])
    >>> ga2equ([359.9443056, -0.0461944444]).round(3)
    array([ 266.417,  -29.008])
    """
    l = np.radians(l)
    b = np.radians(b)

    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = np.radians(192.859508)
    pole_dec = np.radians(27.128336)
    posangle = np.radians(122.932-90.0)
    
    # North galactic pole (B1950)
    #pole_ra = radians(192.25)
    #pole_dec = radians(27.4)
    #posangle = radians(123.0-90.0)
    
    ra = np.arctan2( (np.cos(b)*np.cos(l-posangle)), (np.sin(b)*np.cos(pole_dec) - np.cos(b)*np.sin(pole_dec)*np.sin(l-posangle)) ) + pole_ra
    dec = np.arcsin( np.cos(b)*np.cos(pole_dec)*np.sin(l-posangle) + np.sin(b)*np.sin(pole_dec) )
    
    return np.degrees(ra), np.degrees(dec)

from astropy import units as u
from astropy.coordinates import SkyCoord

def ga2equ_old(l,b):
    c=SkyCoord(l=l*u.degree,b=b*u.degree,frame='galactic')
    return c.fk5.ra.degree, c.fk5.dec.degree

def xy_funk(l,b):
    l=np.where(l>180,l-360,l)    
    lam=-np.radians(l)
    phi=np.radians(b)

    theta=phi
    for i in range(10):
        theta=theta-(2*theta+np.sin(2*theta)-3.1415*np.sin(phi))/(2+2*np.cos(2*theta))

    return 2.0/3.1415*lam*np.cos(theta), np.sin(theta)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(90.-decl),np.radians(RA))

rah,ram,ras,ded,dem,des,l, b, v = np.loadtxt("hzoa.txt",usecols=[1,2,3,4,5,6,7,8,10],unpack=True,skiprows=65)
ra2=(rah+(ram+ras/60)/60)*15
dec2=np.where(ded<0,ded-(dem+des/60)/60,ded+(dem+des/60)/60)
ra,dec=ga2equ(l,b)
np.savetxt('testconv.dat',np.transpose([ra,ra2,dec,dec2]))
ii_hzoa=DeclRaToIndex(dec,ra,256)
x,y =xy_funk(l,b)

map=hp.read_map("2MPZ.gz_0.01_0.02_smoothed_inf.fits.gz",0)
print(len(l))

mask=hp.read_map("prod_mask.fits.gz",0)
mask=hp.pixelfunc.ud_grade(mask,nside_out = 256, order_in = 'RING', order_out = 'RING')

mapf=map*(1-mask)

print('mean on HZOA %g ; mean in mask %g' %
      (np.mean(mapf[ii_hzoa]),np.mean(mapf)))
    
sig=hp.read_map("text/sigmamap.fits.gz",0)
wmapf=mapf/sig**2
wf=(1-mask)/sig**2
print('sigma mean on HZOA %g ; mean in mask %g' %
      (np.sum(wmapf[ii_hzoa])/np.sum(wf[ii_hzoa]),
       np.sum(wmapf)/np.sum(wf)))

nn=len(l)
v=[]
w=[]

for i in range(10000):
    longr=np.random.uniform(low=212-360,high=36,size=nn)
    latr=np.random.uniform(low=-5,high=5,size=nn)
    rar,decr=ga2equ(longr,latr)
    ii_rand=DeclRaToIndex(decr,rar,256)
    v.append(np.mean(mapf[ii_rand]))
    w.append(np.sum(wmapf[ii_rand])/np.sum(wf[ii_rand]))
    
np.savetxt('dist_hzoa.txt',np.transpose([np.sort(v),(np.arange(len(v))+1.0)/len(v)]))
np.savetxt('sdist_hzoa.txt',np.transpose([np.sort(w),(np.arange(len(w))+1.0)/len(w)]))

print('v: %g %g' % (np.mean(v),np.std(v)))
print('w: %g %g' % (np.mean(w),np.std(w)))


ra2,dec2,l2, b2 = np.loadtxt("leda.txt",usecols=[3,4,7,8],unpack=True,skiprows=1)
x2,y2 =xy_funk(l2,b2)
# ra2,dec2=ga2equ([l2,b2])
ii_leda=DeclRaToIndex(dec2,ra2,256)
print('mean on LEDA %g' %
      np.mean((map*(1-mask))[ii_leda]))

hp.mollview(np.where(map>0.01,np.log10(map),-2)+2,coord=['C','G'],title='')
hp.graticule()
# keep=np.logical_and(v>4500,v<5500)
# x=x[keep]
# y=y[keep]
plt.scatter(x,y,color='magenta',marker='*')
# plt.scatter(x2,y2,c='b',s=0.2)
plt.show()
