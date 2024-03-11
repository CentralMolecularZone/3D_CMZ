import pylab as pl
import numpy as np
from numpy import pi, sqrt, cos, sin, arctan2, linspace, log, exp, sinh, arcsinh, arccos, arcsin

##########################
# global variables
##########################

sgr_a_coords     = np.array([359.94423568-360,-0.04616002]) # [deg]
sun_UVW_motion   = np.array([10.0,11.0,7.0]) # [km/s] (U,V,W) Solar motion with respect to  the  Local  Standard  of  Rest. See Sect 5.3.3 of Bland-Hawthorn & Gerhard 2016 (https://arxiv.org/pdf/1602.07702.pdf)
sun_vcirc        = np.array([238.0]) # [km/s] circular rotatio velocity at the Sun 
Rsun             = 8.2   # [kpc] distance of Sun from the Galactic Centre 
zsun             = 0.025 # [kpc] distance of Sun from the Galactic Plane

sun_total_motion = sun_UVW_motion + np.array([0.0,sun_vcirc[0],0.0])

###################################
# unit conversions
###################################

rad_to_mas       = (60*60*180/pi*1e3)    # converts 1 radians to milliarcseconds
year_to_sec      = (60*60*24*365)        # converts 1 year into seconds
kpc_to_km        = 3.0856775814913673e16 # converts 1 kpc to km

###################################
# xyz are cartesian coordinates centered at the GC
# XYZ are cartesian coordinates centered at the Sun position. X points to the GC, Y points in the -(X \times z) direction, Z points is obtained from the right-hand rule from the previous two.
# lbr are spherical coordinates corresponding to XYZ centered at the Sun position (i.e., the usual Galactic coordinates)
# the Sun position is assumed to be (xsun,ysun,zsun)
# the Sun velocity is assumed to be (vxsun,vysun,vzsun)
# the following functions convert back and forth between the above coordinates
###################################

def XYZhat(xsun,ysun,zsun):
  # returns: Xhat, Yhat, Zhat components in the xyz frame
  rsun = sqrt(xsun**2+ysun**2+zsun**2)
  Rsun = sqrt(xsun**2+ysun**2)
  sintheta = Rsun/rsun
  costheta = zsun/rsun
  sinphi   = ysun/Rsun
  cosphi   = xsun/Rsun 
  Xhat = np.array([-cosphi*sintheta,-sinphi*sintheta,-costheta])
  Yhat = np.array([+sinphi,         -cosphi,         +0.0     ])
  Zhat = np.array([-cosphi*costheta,-sinphi*costheta,+sintheta])
  return Xhat, Yhat, Zhat

def rblhat(l,b,r):
  # returns: rhat, bhat, lhat components in the xyz frame
  theta = pi/2 - b
  rhat  = [+sin(theta)*cos(l),+sin(theta)*sin(l),+cos(theta)]
  bhat  = [-cos(theta)*cos(l),-cos(theta)*sin(l),+sin(theta)] 
  lhat  = [-sin(l),cos(l),0]
  return rhat, bhat, lhat

def xyz2XYZ(x,y,z,vx,vy,vz,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
  Xhat,Yhat,Zhat = XYZhat(xsun,ysun,zsun)
  Deltax  = x-xsun
  Deltay  = y-ysun
  Deltaz  = z-zsun
  Deltavx = vx-vxsun
  Deltavy = vy-vysun
  Deltavz = vz-vzsun
  X  = Deltax *Xhat[0] + Deltay *Xhat[1] + Deltaz *Xhat[2]
  Y  = Deltax *Yhat[0] + Deltay *Yhat[1] + Deltaz *Yhat[2]
  Z  = Deltax *Zhat[0] + Deltay *Zhat[1] + Deltaz *Zhat[2]
  vX = Deltavx*Xhat[0] + Deltavy*Xhat[1] + Deltavz*Xhat[2]
  vY = Deltavx*Yhat[0] + Deltavy*Yhat[1] + Deltavz*Yhat[2]
  vZ = Deltavx*Zhat[0] + Deltavy*Zhat[1] + Deltavz*Zhat[2]
  return X,Y,Z,vX,vY,vZ

def XYZ2xyz(X,Y,Z,vX,vY,vZ,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
  Xhat,Yhat,Zhat = XYZhat(xsun,ysun,zsun)
  Deltax  = X *Xhat[0] + Y *Yhat[0] + Z *Zhat[0]
  Deltay  = X *Xhat[1] + Y *Yhat[1] + Z *Zhat[1]
  Deltaz  = X *Xhat[2] + Y *Yhat[2] + Z *Zhat[2]
  Deltavx = vX*Xhat[0] + vY*Yhat[0] + vZ*Zhat[0]
  Deltavy = vX*Xhat[1] + vY*Yhat[1] + vZ*Zhat[1]
  Deltavz = vX*Xhat[2] + vY*Yhat[2] + vZ*Zhat[2]
  x  = Deltax + xsun
  y  = Deltay + ysun
  z  = Deltaz + zsun
  vx = Deltavx + vxsun
  vy = Deltavy + vysun
  vz = Deltavz + vzsun
  return x,y,z,vx,vy,vz

def XYZ2lbr(X,Y,Z,vX,vY,vZ):
  r     = sqrt(X**2+Y**2+Z**2)
  l     = arctan2(Y,X)
  theta = arccos(Z/r)
  b     = pi/2 - theta
  rhat, bhat, lhat = rblhat(l,b,r)
  vr = vX*rhat[0] + vY*rhat[1] + vZ*rhat[2]
  vb = vX*bhat[0] + vY*bhat[1] + vZ*bhat[2]
  vl = vX*lhat[0] + vY*lhat[1] + vZ*lhat[2]
  return l,b,r,vl,vb,vr

def lbr2XYZ(l,b,r,vl,vb,vr):
  theta = pi/2 - b
  rhat, bhat, lhat = rblhat(l,b,r)
  X  = r*sin(theta)*cos(l)
  Y  = r*sin(theta)*sin(l)
  Z  = r*cos(theta)
  vX = vr*rhat[0] + vl*lhat[0] + vb*bhat[0]
  vY = vr*rhat[1] + vl*lhat[1] + vb*bhat[1]
  vZ = vr*rhat[2] + vl*lhat[2] + vb*bhat[2]
  return X,Y,Z,vX,vY,vZ

def xyz2lbr(x,y,z,vx,vy,vz,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
  X,Y,Z,vX,vY,vZ = xyz2XYZ(x,y,z,vx,vy,vz,xsun,ysun,zsun,vxsun,vysun,vzsun)
  l,b,r,vl,vb,vr = XYZ2lbr(X,Y,Z,vX,vY,vZ)
  return l,b,r,vl,vb,vr

def lbr2xyz(l,b,r,vl,vb,vr,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
  X,Y,Z,vX,vY,vZ = lbr2XYZ(l,b,r,vl,vb,vr)
  x,y,z,vx,vy,vz = XYZ2xyz(X,Y,Z,vX,vY,vZ,xsun,ysun,zsun,vxsun,vysun,vzsun)
  return x,y,z,vx,vy,vz

####################
# proper motions
#####################

def vlb_2_mulb(r,v_l,v_b):
    # assumes that velocities are in km/s and proper motions in mas/yr
    mu_l = (v_l/r)*(rad_to_mas*year_to_sec/kpc_to_km)
    mu_b = (v_b/r)*(rad_to_mas*year_to_sec/kpc_to_km)
    return mu_l, mu_b

def mulb_2_vlb(r,mu_l,mu_b):
    # assumes that velocities are in km/s and proper motions in mas/yr
    v_l = (mu_l*r)/(rad_to_mas*year_to_sec/kpc_to_km)
    v_b = (mu_b*r)/(rad_to_mas*year_to_sec/kpc_to_km)
    return v_l, v_b

def subtract_sun_motion(l,b,r,vl,vb,vr,vXsun=0.0,vYsun=2.2,vZsun=0.0):
  # take (l,b,v,vl,vb,vr) coordinates and subtract the Sun motion 
  # i.e., update the velocities to the values that you would observe at the Sun if the Sun were not moving with respect to the GC
  # note vXsun,vYsun and vZsun are in the XYZ coordinate system!
  X,Y,Z,vX,vY,vZ = lbr2XYZ(l,b,r,vl,vb,vr)
  vX = vX + vXsun
  vY = vY + vYsun
  vZ = vZ + vZsun
  l,b,r,vl,vb,vr =  XYZ2lbr(X,Y,Z,vX,vY,vZ)
  return l,b,r,vl,vb,vr 

####################
# sanity checks
#####################

def lbr2xyz_sanity_check(figure=True):
    xsun,  ysun,  zsun  = 0.0,  -8.0, 0.0
    vxsun, vysun, vzsun = -220.0, 0.0, 0.0
    
    # R, theta       = 3.0, np.radians(270)  
    # x,  y,  z      = R*cos(theta), R*sin(theta), 10.0
    # vx, vy, vz     = 0*sin(theta), 0*cos(theta), 10.0

    x,  y,  z      = np.random.rand(3)*10 - 5
    vx, vy, vz     = np.random.rand(3)*200 - 100

    l,b,r,vl,vb,vr       = xyz2lbr(x,y,z,vx,vy,vz,xsun=xsun,ysun=ysun,zsun=zsun,vxsun=vxsun,vysun=vysun,vzsun=vzsun)
    x1,y1,z1,vx1,vy1,vz1 = lbr2xyz(l,b,r,vl,vb,vr,xsun=xsun,ysun=ysun,zsun=zsun,vxsun=vxsun,vysun=vysun,vzsun=vzsun)

    print("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g"%(x,y,z,vx,vy,vz))
    print("l=%g, b=%g, r=%g, vl=%g, vb=%g, vr=%g"%(l,b,r,vl,vb,vr))

    print("original vs check vs difference: %g %g %g "%(x,x1,x-x1))
    print("original vs check vs difference: %g %g %g "%(y,y1,y-y1))
    print("original vs check vs difference: %g %g %g "%(z,z1,z-z1))
    print("original vs check vs difference: %g %g %g "%(vx,vx1,vx-vx1))
    print("original vs check vs difference: %g %g %g "%(vy,vy1,vy-vy1))
    print("original vs check vs difference: %g %g %g "%(vz,vz1,vz-vz1))

    if(figure):
        fig, ax = pl.subplots(figsize=(10,5))
        ax.plot(xsun,ysun,'.',color='b',markersize=10)
        ax.plot(x,y,'.',color='r',markersize=10)
        ax.quiver(x,y,vx,vy,color='r',scale_units='x',scale=50.0)
        ax.quiver(xsun,ysun,vxsun,vysun,color='k',scale_units='x',scale=50.0)
        ax.set_xlabel(r'$x\, {\rm [deg]}$',fontsize=22)
        ax.set_ylabel(r'$y\, {\rm [deg]}$',fontsize=22)
        ax.set_aspect('equal')
        ax.set_xlim(-10,10)
        ax.set_ylim(-10,10)
        ax.grid()
        ax.tick_params(labelsize=18)
        pl.show()

def subtract_sun_motion_sanity_check(figure=True):
    xsun,  ysun,  zsun  = 0.0,  -8.0, 0.0
    vxsun, vysun, vzsun = -220.0, 0.0, 50.0

    dummy,dummy,dummy,vXsun,vYsun,vZsun = xyz2XYZ(0,0,0,0,0,0,xsun=xsun,ysun=ysun,zsun=zsun,vxsun=vxsun,vysun=vysun,vzsun=vzsun)
    vXsun,vYsun,vZsun = -vXsun,-vYsun,-vZsun

    R, theta       = 3.0, np.radians(270)  
    x,  y,  z      = R*cos(theta), R*sin(theta), 0.0
    vx, vy, vz     = 220*sin(theta), 0*cos(theta), 50.0
    
    X,Y,Z,vX,vY,vZ  = xyz2XYZ(x,y,z,vx,vy,vz,xsun=xsun,ysun=ysun,zsun=zsun,vxsun=vxsun,vysun=vysun,vzsun=vzsun)
    # X,Y,Z,vX,vY,vZ = lbr2XYZ(l,b,r,vl,vb,vr)

    l, b, r, vl, vb, vr  = xyz2lbr(x,y,z,vx,vy,vz,xsun=xsun,ysun=ysun,zsun=zsun,vxsun=vxsun,vysun=vysun,vzsun=vzsun)
    l1,b1,r1,vl1,vb1,vr1 = subtract_sun_motion(l,b,r,vl,vb,vr,vXsun=vXsun,vYsun=vYsun,vZsun=vZsun)
  
    print("vxsun=%g, vysun=%g, vzsun=%g"%(vxsun,vysun,vzsun))
    print("vXsun=%g, vYsun=%g, vZsun=%g"%(vXsun,vYsun,vZsun))

    print("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g"%(x,y,z,vx,vy,vz))
    print("X=%g, Y=%g, Z=%g, vX=%g, vY=%g, vZ=%g"%(X,Y,Z,vX,vY,vZ))
    print("l=%g, b=%g, r=%g, vl=%g, vb=%g, vr=%g"%(l,b,r,vl,vb,vr))
    print("l1=%g, b1=%g, r1=%g, vl1=%g, vb1=%g, vr1=%g"%(l1,b1,r1,vl1,vb1,vr1))

    if(figure):
        fig, ax = pl.subplots(figsize=(10,5))
        ax.plot(xsun,ysun,'.',color='b',markersize=10)
        ax.plot(x,y,'.',color='r',markersize=10)
        ax.quiver(x,y,vx,vy,color='r',scale_units='x',scale=50.0)
        ax.quiver(xsun,ysun,vxsun,vysun,color='k',scale_units='x',scale=50.0)
        ax.set_xlabel(r'$x\, {\rm [deg]}$',fontsize=22)
        ax.set_ylabel(r'$y\, {\rm [deg]}$',fontsize=22)
        ax.set_aspect('equal')
        ax.set_xlim(-10,10)
        ax.set_ylim(-10,10)
        ax.grid()
        ax.tick_params(labelsize=18)
        pl.show()


