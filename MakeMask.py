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
            'A FITS file to output the mask map in HEALPIX format'
        ),
    )
    arguments = vars(parser.parse_args())
    return arguments

def MakeMask(galmapfile,maskfile)
    randmap=hp.read_map(galmapfile,0)
    hp.write_map(inffile, yt)


#------------------------------------------------------------------------------
# main
#
def _main():
    """
    This is the main routine.
    """

    args=_parse_command_line_arguments()
    
    MakeMask(args['sky-map'],args['mask-map'])

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    _main()



