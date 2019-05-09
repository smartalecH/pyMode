from pyMode.materials import Si, SiO2
import pyMode as pm
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import pscan

m = pm.PML(pm.Location.N)
e = pm.PML(pm.Location.E)
s = pm.Magnetic(pm.Location.S)

#makes grid that gives locations for grid for FDTD
#takes in set locations, and densities at those locations
def makeGrid(xs,hs):
    grid = [min(xs)]
    interp = interp1d(xs, hs, kind='linear')
    while grid[-1] < max(xs):
        h = interp(grid[-1])
        grid.append(grid[-1] + h)
    return np.array(grid)

#boundaries conditions
boundaries = []#[m,e,f]

width = 0.5        #actual width of (top) of waveguide
thickness = 0.22   #actual thickness (height) of waveguide
xWidth    = 2      #width of full testing location
yWidth    = 1      #height of full testing location
dcladding = 51.1e-3 #line density to use for cladding
dcore     = [21.1e-3,21.1e-3,10.1e-3,5.1e-3, 2.1e-3] #density to use in core
k0        = np.zeros((len(dcore),)) #save answers here
for k in range(len(dcore)):

    #it's a trapezoid, so calculate numbers
    sidewallAngle = 85
    sidewallAngle_radians = sidewallAngle / 180 * np.pi
    bottomFace = width + 2*(thickness/np.tan(sidewallAngle_radians)) #size of the bottom of waveguid
    print(bottomFace)

    #make waveguide centered at (0,0)
    waveguide1 = pm.Trapezoid(
        center=pm.Vector3(0,0),
        topFace = width,
        thickness = thickness,
        sidewallAngle = sidewallAngle_radians,
        core = Si,
        cladding=SiO2,
        rc = 0.0
    )

    xLocs = [0,bottomFace/2,bottomFace/2+0.1,xWidth/2] #where the densities should change
    xVals = [dcore[k],dcore[k],dcladding,dcladding] #densities at those locations (should be 1 more than above)
    xLeft = makeGrid(xLocs,xVals) #makes the actual grid from left side
    xx    = xLeft - xLeft[-1]
    xx    = np.concatenate((-(np.flip(xLeft[1:-1],0)),xLeft)) #use if not using boundary condition
    xx = xLeft

    #same thing as before
    yLocs = [0,thickness/2,thickness/2+0.1,yWidth/2]
    yVals = [dcore[k],dcore[k],dcladding,dcladding]
    yLeft = makeGrid(yLocs,yVals)
    yy    = yLeft - yLeft[-1]
    yy    = np.concatenate((-(np.flip(yLeft[1:-1],0)),yLeft))
    #yy    = yLeft


    geometry = [waveguide1] #include more values here if there's more than 1 geometry
    wavelength = 1.5
    numModes = 1
    radius = 0

    #start simulation
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

    #run it
    sim.run()

    #get answers
    eps = sim.getEps()
    k0[k], Ex = sim.getFieldComponent('hr',0)


plt.figure()
plt.semilogx(dcore,k0/(2*np.pi/wavelength),'o')
plt.savefig('results.png')





