import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import math as mt
import sys

def IndexToDeclRa(NSIDE,index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(mt.pi/2.0-theta),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(90.-decl),np.radians(RA))

for d in sys.argv[1:]:
    try:
        output_map = hp.read_map("%s/final_map.fits.gz" % d,0)
        input_map = hp.read_map("%s/rand_map.fits.gz" % d,0)
        nside = hp.get_nside(input_map)
        comb_mask = 1-hp.read_map("%s/mask.fits.gz" % d,0)

        ndata=len(comb_mask[comb_mask==1])
        # apply the mask
        input_map=input_map*comb_mask
        
        input_map=(input_map-np.sum(input_map)/ndata)*comb_mask
        norm_input=np.sum(input_map*input_map)
        
        output_map=output_map*comb_mask
        output_map=(output_map-np.sum(output_map)/ndata)*comb_mask
        norm_output= np.sum(output_map*output_map)
        
        print('%d %g %g %s' %
              (ndata,np.sum(input_map*input_map)/(norm_input*norm_input)**0.5,
               np.sum(output_map*input_map)/(norm_input*norm_output)**0.5,d))
    except:
        pass
