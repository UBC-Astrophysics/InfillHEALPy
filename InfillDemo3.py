#!/usr/bin/env python
#
# InfileDemo3.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script creates a simple mask, calculates the alm matrix, alm matrix inverse, saves them
# and demonstrates on a test image
#
# usage: python InfillDemo.py
#
#
# Questions: heyl@phas.ubc.ca
#
#     Copyright 2016, Elisa Antolini and Jeremy Heyl
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from argparse import ArgumentParser

import math as mt
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pyfits

def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(mt.pi/2.0-theta),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(90.-decl),np.radians(RA))


nside=64
mask=np.zeros(hp.nside2npix(nside))
i=np.arange(hp.nside2npix(nside))
decl,ra=IndexToDeclRa(nside,i)
mask=abs(decl)+(180-abs(ra-180))/12-10
mask2=np.where(mask>0,np.exp(-0.1/(mask*mask)),0)
mask=mask2
#mask=abs(decl)+abs(ra)/15>12
#mask=abs(decl)>10
if True:
    #Create the pixel map
    lmax=32
    mmax=lmax
    almsize=mmax*(2*lmax+1-mmax)/2+lmax+1
    almmatrix=[]
    for i in range(almsize):
        x=np.zeros(almsize)
        x[i]=1
        alm=x.astype(np.complex128)
        map=hp.sphtfunc.alm2map(alm,nside)
        
        if False:
            hp.mollview(map,coord='C',rot = [0,0.3],
                        title='Test Plot', unit='prob', xsize=1024)
            hp.graticule()
            plt.show()
            
        maskedmap=mask*map
        almnew=hp.sphtfunc.map2alm(maskedmap,lmax=lmax)
        almmatrix.append(almnew)

        # the pseudo inverse is better
    almpinvmatrix=np.linalg.pinv(np.transpose(almmatrix))
    np.savez_compressed("maskdata.npz",mask=mask,almmatrix=almmatrix,lmax=lmax,mmax=mmax,almpinvmatrix=almpinvmatrix)
else:
    pass


#Create random map
if True:
    lmax2=16
    mmax2=lmax2
    almsize2=mmax2*(2*lmax2+1-mmax2)/2+lmax2+1
    realrandalm=np.random.uniform(size=almsize2)
    imagrandalm=np.random.uniform(size=almsize2)
    randmap=hp.sphtfunc.alm2map(realrandalm+1j*imagrandalm,nside)
else:
#load galaxy map
    randmap=hp.read_map("2MPZ.gz_0.031_0.051_smoothed.fits",0)
    randmap=hp.pixelfunc.ud_grade(randmap,nside_out = nside, order_in = 'RING', order_out = 'RING')

hp.mollview(mask,coord='C',rot = [0,0.3],
            title='Mask', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("Mask.png")
plt.show()

hp.mollview(randmap,coord='C',rot = [0,0.3],
            title='RandMap', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("RandMap.png")
plt.show()

maskedrandmap=mask*randmap
hp.mollview(maskedrandmap,coord='C',rot = [0,0.3],
            title='MaskedRandMap', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("MaskedRandMap.png")
plt.show()

almmaskedrandmap=hp.sphtfunc.map2alm(maskedrandmap,lmax=lmax,mmax=mmax)
print(len(almmaskedrandmap))

hp.mollview(hp.sphtfunc.alm2map(almmaskedrandmap,nside),coord='C',rot = [0,0.3],
            title='LowPassMaskedRandMap', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("LowPassMaskedRandMap.png")
plt.show()

lmax2=lmax
if lmax2<lmax:
    mmax2=lmax2
    almsize2=mmax2*(2*lmax2+1-mmax2)/2+lmax2+1
    # truncatedalmmatrix=np.matrix(almmatrix[0:almsize2])
    # almvector=np.matrix(almmaskedrandmap)
    # almunmaskedrandmap=np.linalg.pinv(truncatedalmmatrix)*almvector
    almunmaskedrandmap=np.dot(np.linalg.pinv(np.transpose(almmatrix[0:almsize2])),almmaskedrandmap)
else:
    # use the pe
    # almunmaskedrandmap=np.concatenate((np.dot(almmaskedrandmap[0:len(alminvmatrix)],alminvmatrix),almmaskedrandmap[len(alminvmatrix):]))
    #almunmaskedrandmap=np.dot(almmaskedrandmap[0:len(alminvmatrix)],alminvmatrix)
    almunmaskedrandmap=np.dot(almpinvmatrix,almmaskedrandmap)
    #almunmaskedrandmap=np.dot(almmatrix,almmaskedrandmap[0:len(almmatrix)])
    #almunmaskedrandmap=almmaskedrandmap[0:len(almmatrix)]

rawunmaskedrandmap=hp.sphtfunc.alm2map(almunmaskedrandmap,nside)

hp.mollview(rawunmaskedrandmap,coord='C',rot = [0,0.3],
            title='RawUnMaskedRandMap', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("RawUnMaskedRandMap.png")
plt.show()

unmaskedrandmap=mask*maskedrandmap+(1-mask)*rawunmaskedrandmap
hp.mollview(unmaskedrandmap,coord='C',rot = [0,0.3],
            title='UnMaskedRandMap', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("UnMaskedRandMap.png")
plt.show()

diffmap=unmaskedrandmap-randmap
hp.mollview(diffmap,coord='C',rot = [0,0.3],
            title='DiffMap', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("DiffMap.png")
plt.show()




