#!/usr/bin/env python
#
# CombineMap.py
#
# Elisa Antolini
# Jeremy Heyl
# UBC Southern Observatory
#
# This script loads a simple mask, image and infilled image, and rescales the
# infilled region so that the mean and dispersion of the infilled region equals
# that of the rest of the map.
#
# usage: CombineMap.py [-h] sky-map mask-map infilled-map combined-map
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


def _parse_command_line_arguments():
    """
    Parse and return command line arguments
    """
    parser = ArgumentParser(
        description=(
            'Command-line tool to Infill a galaxy map within a masked region'
        ),
    )
    parser.add_argument(
        'sky-map',
        type=str,
        help=(
            'A FITS file containing the Galaxy map in HEALPIX format'
        ),
    )
    parser.add_argument(
        'mask-map',
        type=str,
        help=(
            'A FITS file containing the mask map in HEALPIX format'
        ),
    )
    parser.add_argument(
        'infilled-map',
        type=str,
        help=(
            'A FITS file for the output infilled map'
        ),
    )
    parser.add_argument(
        'combined-map',
        type=str,
        help=(
            'A FITS file for the output combined map'
        ),
    )
    arguments = vars(parser.parse_args())
    return arguments

def CombineMap(galmapfile,maskfile,inffile,combfile):
    randmap=hp.read_map(galmapfile,0)
    mask=hp.read_map(maskfile,0)

    nside=hp.get_nside(randmap)

    mask=hp.pixelfunc.ud_grade(mask,nside_out = nside, order_in = 'RING', order_out = 'RING')
    
    infmap=hp.read_map(inffile,0)

    # calculate the mean and sigma of the unmasked galmap file
    maskedrandmap=mask*randmap
    meanmap=np.mean(maskedrandmap)/np.mean(mask)
    sigmamap=np.std(maskedrandmap)/mt.sqrt(np.mean(mask*mask))

    # calculate the mean and sigma of the masked region in the infill map
    onemmask=1-mask
    maskedinfmap=onemmask*infmap
    meaninfmap=np.mean(maskedinfmap)/np.mean(onemmask)
    sigmainfmap=np.std(maskedinfmap)/mt.sqrt(np.mean(onemmask**2))

    combinedmap=maskedrandmap+onemmask*((maskedinfmap-meaninfmap)/sigmainfmap+meanmap)*sigmamap
    hp.write_map(combfile,combinedmap)

#------------------------------------------------------------------------------
# main
#
def _main():
    """
    This is the main routine.
    """

    args=_parse_command_line_arguments()
    
    CombineMap(args['sky-map'],args['mask-map'],args['infilled-map'],
           args['combined-map'])

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    _main()



