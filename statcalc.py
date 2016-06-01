import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
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

if False:
    hp.mollview(input_map*comb_mask,coord=['C'],title='',min=0,max=0.65)
    plt.show()
    hp.mollview(output_map*comb_mask,coord=['C'],title='',min=0,max=0.65)
    plt.show()

# calculate the cl of the input_map to generate new ones
if False:
    cl=hp.anafast(input_map)

ndata=len(comb_mask[comb_mask==1])
print("ndata = %d" % ndata)
# apply the mask
input_map=input_map*comb_mask

input_map=(input_map-np.sum(input_map)/ndata)*comb_mask
norm_input=np.sum(input_map*input_map)

output_map=output_map*comb_mask
output_map=(output_map-np.sum(output_map)/ndata)*comb_mask
norm_output= np.sum(output_map*output_map)

print('%d %g' %
      (-2,np.sum(input_map*input_map)/(norm_input*norm_input)**0.5))
print('%d %g' %
      (-1,np.sum(output_map*input_map)/(norm_input*norm_output)**0.5))

if True:
    index=np.arange(len(input_map))
    theta,phi= IndexToDeclRa(nside,index)
    for ang in np.arange(3600)*0.1:
        index_new=DeclRaToIndex(theta,phi+ang,nside)
        comb_mask_new=comb_mask[index_new]*comb_mask
        ndata_new=len(comb_mask_new[comb_mask_new==1])
        output_map_new=output_map[index_new]*comb_mask_new
        output_map_new=(output_map_new-np.sum(output_map_new)/ndata_new)*comb_mask_new
        input_map_new=input_map*comb_mask_new
        input_map_new=(input_map_new-np.sum(input_map_new)/ndata_new)*comb_mask_new
        print('%g %d %g' %
              (ang,ndata_new,np.sum(output_map_new*input_map_new)/(
                  np.sum(input_map_new*input_map_new)*
                  np.sum(output_map_new*output_map_new))**0.5))
    
    

if False:
    for cc in range(1000):
        newmap=hp.synfast(cl,nside=nside)
        newmap=newmap*comb_mask 
        newmap=(newmap-np.sum(newmap)/ndata)*comb_mask
        print('%d %g' %
              (cc,
               np.sum(output_map*newmap)/
               (np.sum(newmap*newmap)*norm_output)**0.5))

if False:
    itest=input_map[comb_mask==1]
    otest=output_map[comb_mask==1]
    print('%d %g' %
          (-1,np.sum(otest*itest)/(norm_input*norm_output)**0.5))
    for cc in range(1000):
        np.random.shuffle(otest)
        print('%d %g' %
          (cc,np.sum(otest*itest)/(norm_input*norm_output)**0.5))
