import healpy as hp
import sys
import matplotlib.pyplot as plt
import numpy 
for rr in sys.argv[1:]:
    map=hp.read_map(rr,0)
    hp.mollview(map,coord=['C','G'],rot = [0,0.3],title='',norm='hist'),
    hp.graticule()
    plt.show()
    plt.close()
