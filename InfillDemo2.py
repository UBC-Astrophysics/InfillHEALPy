#!/usr/bin/env python
#
# InfillDemo2.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script loads a mask, the alm matrix and alm matrix inverse from a file 
# and demonstrates on a test image
#
# usage: python InfillDemo2.py
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

def _CreateRandomMap(lmax2,nside):
    #Create random map
    mmax2=lmax2
    almsize2=mmax2*(2*lmax2+1-mmax2)/2+lmax2+1
    realrandalm=np.random.uniform(size=almsize2)
    imagrandalm=np.random.uniform(size=almsize2)

    return hp.sphtfunc.alm2map(realrandalm+1j*imagrandalm,nside)

def _doInfillDemo(randmap,mask,alminvmatrix,nside,lmax,mmax,lmax2=16):

    # 1. mask the image
    maskedrandmap=mask*randmap

    # 2. calculate the alm of the masked image
    almmaskedrandmap=hp.sphtfunc.map2alm(maskedrandmap,lmax=lmax,mmax=mmax)
    
    # 3. calculate the best guess of the alm of the orignal alm using the alminvmatrix
    almunmaskedrandmap=np.dot(alminvmatrix,almmaskedrandmap)

    # 4. calculate the ifft of the unmasked alm
    rawunmaskedrandmap=hp.sphtfunc.alm2map(almunmaskedrandmap,nside)

    # 5. only use the ifft in the masked area to create the unmasked map
    unmaskedrandmap=mask*maskedrandmap+(1-mask)*rawunmaskedrandmap
    
    # now do the plots
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

    hp.mollview(maskedrandmap,coord='C',rot = [0,0.3],
                title='MaskedRandMap', unit='prob', xsize=1024)
    hp.graticule()
    plt.savefig("MaskedRandMap.png")
    plt.show()

    hp.mollview(hp.sphtfunc.alm2map(almmaskedrandmap,nside),coord='C',rot = [0,0.3],
                title='LowPassMaskedRandMap', unit='prob', xsize=1024)
    hp.graticule()
    plt.savefig("LowPassMaskedRandMap.png")
    plt.show()

    hp.mollview(rawunmaskedrandmap,coord='C',rot = [0,0.3],
                title='RawUnMaskedRandMap', unit='prob', xsize=1024)
    hp.graticule()
    plt.savefig("RawUnMaskedRandMap.png")
    plt.show()

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

def _main():
    a=np.load("maskdata.npz")
    mask=a['mask']
    # calculate nside from the healpix array
    nside=int(np.sqrt(len(mask)/12))
    alminvmatrix=a['alminvmatrix']
    lmax=a['lmax']
    mmax=a['mmax']
    randmap=_CreateRandomMap(16,nside)
    _doInfillDemo(randmap,mask,alminvmatrix,nside,lmax,mmax)
    
#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    _main()



