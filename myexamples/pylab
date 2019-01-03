
import numpy as np

angfac = 180.0/np.pi

# some useful subroutines

# get an angle between [0,2pi]
def mod_two_pi(x):
    twopi = 2.0*np.pi
    y=x
    while (y > twopi):
        y -= twopi;
    while (y < 0.0):
        y += twopi;
    return y

def mod_two_pi_arr(x):
    nvec = np.size(x)
    mvec = x*0.0
    for i in range(0,nvec):
        mvec[i] = mod_two_pi(x[i])
    return mvec

# length of a vector
def len_vec(x,y,z):
    r= np.sqrt(x*x + y*y + z*z)
    return r

# normalize a vector
def normalize_vec(x,y,z):
    r = len_vec(x,y,z)
    return x/r, y/r, z/r

def dotprod(ax,ay,az,bx,by,bz):
    z = ax*bx + ay*by + az*bz  # dot product
    return z

# return cross product of two vectors
def crossprod(ax,ay,az,bx,by,bz):
    cx = ay*bz-az*by;
    cy = az*bx-ax*bz;
    cz = ax*by-ay*bx;
    return cx,cy,cz

# return normalized cross product of two vectors
def crossprod_unit(ax,ay,az,bx,by,bz):
    cx,cy,cz = crossprod(ax,ay,az,bx,by,bz)
    cc = len_vec(cx,cy,cz)
    return cx/cc, cy/cc, cz/cc

# return the vector part of a that is perpendicular to b direction
def aperp(ax,ay,az,bx,by,bz):
    z = ax*bx + ay*by + az*bz  # dot product = ab cos theta
    bmag = len_vec(bx,by,bz)
    #theta = np.acos(z/(amag*bmag))
    cx = ax - z*bx/bmag
    cy = ay - z*by/bmag
    cz = az - z*bz/bmag
    return cx,cy,cz

# return the vector part of a that is parallel to b direction
def apar(ax,ay,az,bx,by,bz):
    z = ax*bx + ay*by + az*bz  # dot product = ab cos theta
    bmag = len_vec(bx,by,bz)
    #theta = np.acos(z/(amag*bmag))
    cx = z*bx/bmag
    cy = z*by/bmag
    cz = z*bz/bmag
    return cx,cy,cz
    

