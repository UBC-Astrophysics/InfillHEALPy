import numpy
import math

def eq2gal(ra,dec):

    rmat = numpy.array([
     [-0.054875539726,-0.873437108010,-0.483834985808],
     [+0.494109453312,-0.444829589425,+0.746982251810],
     [-0.867666135858,-0.198076386122,+0.455983795705]],dtype='d')
    cosb = math.cos(dec)
    v1 = numpy.array([math.cos(ra)*cosb,math.sin(ra)*cosb,math.sin(dec)])
    v2=numpy.dot(rmat,v1)
    x = v2[0]
    y = v2[1]
    z = v2[2]
    r=math.sqrt(x*x+y*y)
    if r == 0.0:
        l = 0.0
    else:
        l = math.atan2(y,x)
    if z == 0.0:
        b = 0.0
    else:
        b = math.atan2(z,r)
    ll = math.fmod(l,2.*math.pi)
    if ll < 0.0: ll = ll + 2.*math.pi

    bb = math.fmod(b,2*math.pi)
    if abs(bb) >= math.pi: print "Ugh!" 
    
    return ll, bb
