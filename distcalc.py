import healpy as hp
import numpy as np
import math as mt

def IndexToDeclRa(NSIDE,index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(mt.pi/2.0-theta),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(90.-decl),np.radians(RA))

output_map = hp.read_map("infill_test.fits.gz",0)
input_map = hp.read_map("2MPZ.gz_0.01_0.1_smoothed.fits.gz",0)
nside = hp.get_nside(input_map)
test_mask = hp.read_map("test_mask.fits.gz",0)
gal_mask = hp.read_map("prod_mask.fits.gz",0)
comb_mask = (1-test_mask)*gal_mask

ndata=len(comb_mask[comb_mask==1])
print("ndata = %d" % ndata)
# apply the mask

imap=input_map[comb_mask>0]
imap=np.where(imap<0,0,imap)
omap=output_map[comb_mask>0]
omap=np.where(omap<0,0,omap)

indexomap=np.argsort(-omap)
ngalaxies=np.cumsum(imap[indexomap])

indexomap=np.argsort(-imap)
ngalperfect=np.cumsum(imap[indexomap])

np.savetxt('test_dist.dat',np.transpose([np.arange(len(imap))*1.0/len(imap),ngalaxies/ngalaxies[-1],ngalperfect/ngalperfect[-1]]))
