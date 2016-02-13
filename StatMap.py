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

