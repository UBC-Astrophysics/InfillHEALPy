#!/usr/bin/env python
#
# PlotMap.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script calculates the mean and dispersion of a series of HEALPix maps
# and outputs the mean and sigma maps
#
# usage: python StatMap.py List_HEALPIX_FITS_Files mean_map sigma_map
#
#
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

import healpy as hp
import numpy as np
import sys

sum=None
cnt=0
for rr in sys.argv[1:-2]:
    map=hp.read_map(rr,0)
    nsidemap=hp.get_nside(map)
    if (sum==None):
        sum=map
        nsidesum=nsidemap
    else:
        if (nsidesum>nsidemap):
            map=hp.pixelfunc.ud_grade(map,nside_out = nsidesum, order_in = 'RING', order_out = 'RING')
        else:
            if (nsidesum<nsidemap):
                sum=hp.pixelfunc.ud_grade(sum,nside_out = nsidemap,
                                          order_in = 'RING',
                                          order_out = 'RING')
                nsidesum=nsidemap
        sum=sum+map
    cnt=cnt+1

mean=sum/cnt

sum=None
for rr in sys.argv[1:-2]:
    map=hp.read_map(rr,0)
    nsidemap=hp.get_nside(map)
    if (sum==None):
        dum=map-mean
        sum=dum*dum
        nsidesum=nsidemap
    else:
        if (nsidesum>nsidemap):
            map=hp.pixelfunc.ud_grade(map,nside_out = nsidesum, order_in = 'RING', order_out = 'RING')
        else:
            if (nsidesum<nsidemap):
                sum=hp.pixelfunc.ud_grade(sum,nside_out = nsidemap,
                                          order_in = 'RING',
                                          order_out = 'RING')
                nsidesum=nsidemap
        dum=map-mean
        sum=sum+dum*dum

sigma=np.sqrt(sum/(cnt-1))

hp.write_map(sys.argv[-2], mean)
hp.write_map(sys.argv[-1], sigma)

