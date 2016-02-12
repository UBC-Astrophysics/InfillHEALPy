#!/usr/bin/env python
#
# Infill.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script creates a simple mask and uses the assumption of sparseness
# in spherical harmonics to infill the mask
#
# usage: python Infill.py
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
import sys


mask=hp.read_map(sys.argv[1],0)

randmap=hp.read_map(sys.argv[2],0)
nside=hp.get_nside(randmap)

mask=hp.pixelfunc.ud_grade(mask,nside_out = nside, order_in = 'RING', order_out = 'RING')

maskedrandmap=mask*randmap

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
    if (i % 5 == 0 and False):
        hp.mollview(yt,coord='C',rot = [0,0.3],
                    title='Reconstruction %d' % i, unit='prob', xsize=1024, norm='hist')
        hp.graticule()
        plt.savefig("recon_%03d.png" % i)
        plt.close()

diffmap=yt-randmap

hp.write_map(sys.argv[3], diffmap)
hp.write_map(sys.argv[4], yt)




