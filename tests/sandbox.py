from pyMode.materials import Si, SiO2
import pyMode as pm
import numpy as np

m = pm.PML(pm.Locations.N)
e = pm.PML(pm.Locations.E)
f = pm.Magnetic(pm.Locations.S)

boundaries = [m,e,f]

waveguide = pm.Rectangle(
    center=pm.Vector3(0,0),
    size = pm.Vector3(1,1),
    core = Si,
    cladding=SiO2
)

geometry = [waveguide]
wavelength = 1.5
numModes = 1
xGrid = np.linspace(-3,3,100)
yGrid = np.linspace(-3,3,100)
radius = 0

sim = pm.Simulation(
    geometry=geometry,
    wavelength=wavelength,
    numModes=numModes,
    xGrid=xGrid,
    yGrid=yGrid,
    radius=radius,
    boundaries = boundaries
    )

sim.run()


