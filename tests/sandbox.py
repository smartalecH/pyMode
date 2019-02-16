from pyMode.materials import Si, SiO2
import pyMode as pm
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

m = pm.PML(pm.Locations.N)
e = pm.PML(pm.Locations.E)
s = pm.Magnetic(pm.Locations.S)

def makeGrid(xs,hs):
    grid = [min(xs)]
    interp = interp1d(xs, hs, kind='linear')
    while grid[-1] < max(xs):
        h = interp(grid[-1])
        grid.append(grid[-1] + h)
    return np.array(grid)


boundaries = []#[m,e,f]

width = 0.5
thickness = 0.22
xWidth    = 2
yWidth    = 1
dcladding = 51.1e-3
dcore     = 1.1e-3


sidewallAngle = 85
sidewallAngle_radians = sidewallAngle / 180 * np.pi
bottomFace = width + 2*(thickness/np.tan(sidewallAngle_radians))
print(bottomFace)


waveguide = pm.Trapezoid(
    center=pm.Vector3(0,0),
    topFace = width,
    thickness = thickness,
    sidewallAngle = sidewallAngle_radians,
    core = Si,
    cladding=SiO2,
    rc = 0.0
)

xLocs = [0,bottomFace/2,bottomFace/2+0.1,xWidth/2]
xVals = [dcore,dcore,dcladding,dcladding]
xLeft = makeGrid(xLocs,xVals)
xx    = xLeft - xLeft[-1]
xx    = np.concatenate((-(np.flip(xLeft[1:-1],0)),xLeft))
xx = xLeft

yLocs = [0,thickness/2,thickness/2+0.1,yWidth/2]
yVals = [dcore,dcore,dcladding,dcladding]
yLeft = makeGrid(yLocs,yVals)
yy    = yLeft - yLeft[-1]
yy    = np.concatenate((-(np.flip(yLeft[1:-1],0)),yLeft))
#yy    = yLeft

plt.figure()
plt.subplot(2,1,1)
plt.plot(xx)
plt.subplot(2,1,2)
plt.plot(yy)
plt.savefig('meshing.png')


geometry = [waveguide]
wavelength = 1.5
numModes = 1
radius = 0

sim = pm.Simulation(
    geometry=geometry,
    wavelength=wavelength,
    numModes=numModes,
    xGrid=xx,
    yGrid=yy,
    radius=radius,
    boundaries = boundaries,
    background = SiO2
    #eigStart = 3
    )

sim.run()

eps = sim.getEps()
k0, Ex = sim.getFieldComponent('er',0)
print(Ex.shape)


X,Y = np.meshgrid(xx,yy)

plt.figure()
plt.contour(X,Y,np.real(Ex))
#plt.pcolor(X,Y,np.real(eps),cmap='gray',alpha=0.2)
plt.tight_layout()
plt.savefig('results.png')



