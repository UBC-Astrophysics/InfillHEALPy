#!/usr/bin/env python
#
# PlotMap.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script plots a bunch of HEALPIX fits file
#
#
# usage: python PlotMap.py List_HEALPIX_FITS_Files
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
import sys
import matplotlib.pyplot as plt
import numpy as np

if (len(sys.argv)<2):
    print("usage: python PlotMap.py List_HEALPIX_FITS_Files")
    exit
    
for rr in sys.argv[1:]:
    map=hp.read_map(rr,0)
#    hp.mollview(np.where(map>0.01,np.log10(map),-2)+2,coord=['C','G'],title='')
    hp.mollview(map,coord=['C'],title='')
    hp.graticule()
    plt.show()
    plt.close()
