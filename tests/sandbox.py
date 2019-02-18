from pyMode.materials import Si, SiO2
import pyMode as pm
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import pscan

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
dcore     = [21.1e-3,21.1e-3,10.1e-3,5.1e-3, 2.1e-3]
k0        = np.zeros((len(dcore),))
for k in range(len(dcore)):

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
    xVals = [dcore[k],dcore[k],dcladding,dcladding]
    xLeft = makeGrid(xLocs,xVals)
    xx    = xLeft - xLeft[-1]
    xx    = np.concatenate((-(np.flip(xLeft[1:-1],0)),xLeft))
    xx = xLeft

    yLocs = [0,thickness/2,thickness/2+0.1,yWidth/2]
    yVals = [dcore[k],dcore[k],dcladding,dcladding]
    yLeft = makeGrid(yLocs,yVals)
    yy    = yLeft - yLeft[-1]
    yy    = np.concatenate((-(np.flip(yLeft[1:-1],0)),yLeft))
    #yy    = yLeft


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
    k0[k], Ex = sim.getFieldComponent('hr',0)


plt.figure()
plt.semilogx(dcore,k0/(2*np.pi/wavelength),'o')
plt.savefig('results.png')





