#!/usr/bin/env python
#
# InfillDemo4.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script creates a simple mask and uses the assumption of sparseness
# in spherical harmonics to infill the mask
#
# usage: python InfillDemo4.py
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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyfits
import eq2gal

def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(mt.pi/2.0-theta),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(90.-decl),np.radians(RA))


nside=64
mask=np.zeros(hp.nside2npix(nside))
i=np.arange(hp.nside2npix(nside))

#Create mask
if True:
    # use RA,Dec coords and mask galaxy
    decl,ra=IndexToDeclRa(nside,i)
    eqfunk=np.vectorize(eq2gal.eq2gal)
    ll,bb=eqfunk(np.radians(ra),np.radians(decl))
    ll=np.degrees(ll)
    bb=np.degrees(bb)
    mask=(abs(bb)+(180-abs(ll-180))/12-10)
#    mask=(abs(bb)-20)
else:
    # use l,b and mask galaxy
    bb, ll=IndexToDeclRa(nside,i)
    mask=(abs(bb)+(180-abs(ll-180))/12-10)
    
#mask=np.where(mask>0,np.exp(-0.1/(mask*mask)),0)
mask=np.where(mask>0,1,0)
hp.write_map("mask.fits.gz", mask)


#Create random map
if False:
    lmax2=16
    mmax2=lmax2
    almsize2=mmax2*(2*lmax2+1-mmax2)/2+lmax2+1
    realrandalm=np.random.uniform(size=almsize2)
    imagrandalm=np.random.uniform(size=almsize2)
    i=np.arange(len(realrandalm))
    weight=1
    randmap=hp.sphtfunc.alm2map(weight*(realrandalm+1j*imagrandalm),nside)
    hp.write_map("rand_map.fits.gz",randmap)
    maskedrandmap=mask*(randmap+np.random.normal(scale=0.2,size=len(randmap)))
    hp.write_map("maskedrand_map.fits.gz",maskedrandmap)
# LOAD RANDOM MAP
if False:
    randmap=hp.read_map("recon_fake/rand_map.fits.gz",0)
    maskedrandmap=hp.read_map("recon_fake/maskedrand_map.fits.gz",0)
# LOAD RANDOM MAP MASK and NOISE
<<<<<<< HEAD
if False:
=======
if True:
>>>>>>> 051706a2d9152003a88d0454cb8dc271be50e22e
    randmap=hp.read_map("../recon_semifake/rand_map.fits.gz",0)
    hp.write_map("rand_map.fits.gz",randmap)
    maskedrandmap=mask*(randmap+np.random.normal(scale=0.02,size=len(randmap)))
    hp.write_map("maskedrand_map.fits.gz",maskedrandmap)
if True:
#load galaxy map
<<<<<<< HEAD
    randmap=hp.read_map("../2MPZ.gz_0.01_0.1_smoothed.fits.gz",0)
=======
    randmap=hp.read_map("2MPZ.gz_0.01_0.1_smoothed.fits.gz",0)
>>>>>>> 051706a2d9152003a88d0454cb8dc271be50e22e
    randmap=hp.pixelfunc.ud_grade(randmap,nside_out = nside, order_in = 'RING', order_out = 'RING')
    meanrandmap=np.mean(randmap)
#    mask[randmap<1e-2*meanrandmap]=0
    maskedrandmap=mask*randmap
# GENERATE RANDOM MAP WITH THE SAME CLS AS LOADED
<<<<<<< HEAD
if True:
=======
if False:
>>>>>>> 051706a2d9152003a88d0454cb8dc271be50e22e
    cl = hp.anafast(randmap)
    randmap=hp.synfast(cl,nside)
    maskedrandmap=mask*randmap
    hp.write_map("rand_map.fits.gz",randmap)
    hp.write_map("maskedrand_map.fits.gz",maskedrandmap)

hp.mollview(mask,coord='C',rot = [0,0.3],
            title='Mask', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("Mask.png")
# plt.show()
plt.close()

hp.mollview(randmap,coord='C',rot = [0,0.3],
            title='Random Map', unit='prob', xsize=1024, norm='hist')
hp.graticule()
plt.savefig("RandMap.png")
# plt.show()
plt.close()


hp.mollview(maskedrandmap,coord='C',rot = [0,0.3],
            title='Masked Random Map', unit='prob', xsize=1024, norm='hist')
hp.graticule()
plt.savefig("MaskedRandMap.png")
# plt.show()
plt.close()



#
# NB: there are more feautures than in the original map
#
lmax2=64
mmax2=64
almsize2=mmax2*(2*lmax2+1-mmax2)/2+lmax2+1

#
# we will threshold and use the strongest features first 
#
# Jerome Bobin, Jean-Luc Starck, Jalal M. Fadili,
# Yassir Moudden, and David L. Donoho (2007)
# Morphological Component Analysis: An Adaptive Thresholding Strategy
# IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 16, NO. 11, 2675
#
# P. Abriala, Y. Mouddena, J.-L. Starck, J. Fadili, J. Delabrouille
# M.K. Nguyen (2008)
# CMB data analysis and sparsity
# Statistical Methodology 5 289-298
#
# Thresholding scheme (Step 1 in Bobin)
#
# 1) Fit using all features
# 2) Sort from strongest to weakest
# 3) Keep only the strongest portion with a fraction=thresh
# 4) thresh increases logarithmically and reaches unity where it
#    iterates further

# set zeroth iteration
i=0

# set initial guess to mean of map (outside of masked region)
yt=0*maskedrandmap+np.sum(maskedrandmap)/np.sum(mask)

for thresh in np.logspace(-3.5,-0.5,200):
    #
    # calculate residual (Step 2a in Bobin et al)
    #
    # we will use just the residual in the unmasked region
    #

    #
    # repeat five times with the same number of coefficients to burn in
    #
    for j in range(1):
        rt=maskedrandmap-yt
    
        #
        # calculate coefficients (Step 2b started in Bobin et al)
        #
        if False:
            alphat=hp.sphtfunc.map2alm(mask*rt+yt,lmax=lmax2,mmax=mmax2)
        else:
            alphat=hp.sphtfunc.map2alm(mask*rt+yt)
            almsize2=len(alphat)
        
        absalphat=np.abs(alphat)
        maxabsalphat=max(absalphat)
        print ('iteration %d %d %g %g %g %d' % (i,j,maxabsalphat,np.sum(absalphat),thresh,int(thresh*almsize2)))
    
        #
        # sort the coefficients from largest to smallest
        #
    
        sortedindex=np.argsort(-absalphat)
    
        #
        # find in the indices of the biggest ones
        #

        bigvalues=sortedindex[0:int(almsize2*thresh)]
        alphatnew=0*alphat
    
        #
        # load the biggest ones into an array with others zeroed out
        #
    
        alphatnew[bigvalues]=alphat[bigvalues]
    
        #
        # apply the thresholds (Step 2b finished)
        #

        alphat=alphatnew
    
        #
        # alphat=np.where(absalphat>thresh,alphat,0)
        # calculate the next guess
        #

        yt=hp.sphtfunc.alm2map(alphat,nside)
    
    i=i+1
#    if (i % 5 == 0):
    if False:
        hp.mollview(yt,coord='C',rot = [0,0.3],
                    title='Reconstruction %d' % i, unit='prob', xsize=1024, norm='hist')
        hp.graticule()
        plt.savefig("recon_%03d.png" % i)
        plt.close()

hp.write_map("final_map.fits.gz",yt)

diffmap=yt-randmap
hp.mollview(diffmap,coord='C',rot = [0,0.3],
            title='DiffMap', unit='prob', xsize=1024)
hp.graticule()
plt.savefig("DiffMap.png")
plt.close()



