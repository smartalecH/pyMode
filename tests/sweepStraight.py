from pyMode.materials import Si, SiO2
import pyMode as pm
import numpy as np
from scipy.interpolate import interp1d
import pscan
from pyMode.materials import Si, SiO2
import argparse
import os
import json
from matplotlib import pyplot as plt

# ---------------------------------------------------------------------------- #
# Parser
# ---------------------------------------------------------------------------- #
def parseInput():
    parser = argparse.ArgumentParser(
        description="Silicon Photonics Waveguide Simulator")
    parser.add_argument("--jobNumber",type=int, required=True)
    parser.add_argument("--numJobs",type=int, required=True)
    parser.add_argument("--dcore",type=float,default=5.1e-3) # resolution in x dimension in microns
    parser.add_argument("--dcladding",type=float,default=30.1e-3) # resolution in y dimension in microns
    parser.add_argument("--xWidth",type=float,default=5)   # domain size in x dimension in microns
    parser.add_argument("--yWidth",type=float,default=3)   # domain size in y dimension in microns
    parser.add_argument("--numModes",type=int,default=1) # number of modes to solve for
    parser.add_argument("--folderName",type=str,default='') # where to output the results to. Note: Must have / at the end.
    parser.add_argument("--filenamePrefix",type=str,default='')   # whether or not this is a testing run for debugging
    
    args = parser.parse_args()
    dict = vars(args)
    for attribute, value in dict.items():
        print('{} : {}'.format(attribute, value))

    # Clean input
    for key,value in dict.items():
        if (dict[key]=="False"):
            dict[key] = False
        elif dict[key]=="True":
            dict[key] = True
        try:
            if dict[key].is_integer():
                dict[key] = int(dict[key])
            else:
                dict[key] = float(dict[key])
        except:
            pass

    kwargs = dict

    # Ensure output directory exists
    folderName = dict['folderName']
    if folderName is not '':
        if not os.path.exists(folderName):
            os.makedirs(folderName)
    # Save dictionary to text file in location specified
    with open(folderName + 'params.txt', 'w') as file:
        file.write(json.dumps(kwargs))
    # Output
    return kwargs
    

# ---------------------------------------------------------------------------- #
# Helper functions
# ---------------------------------------------------------------------------- #

def makeGrid(xs,hs):
    grid = [min(xs)]
    interp = interp1d(xs, hs, kind='linear')
    while grid[-1] < max(xs):
        h = interp(grid[-1])
        grid.append(grid[-1] + h)
    return np.array(grid)

# ---------------------------------------------------------------------------- #
# Sweeping routine
# ---------------------------------------------------------------------------- #

def runSimulation(wavelength,width,thickness,sidewallAngle,dcore,dcladding,xWidth,yWidth,numModes=1):
    
    # Set up the waveguide geometry
    sidewallAngle_radians = sidewallAngle / 180 * np.pi
    bottomFace = width + 2*(thickness/np.tan(sidewallAngle_radians))

    waveguide = pm.Trapezoid(
        center=pm.Vector3(0,0),
        topFace = width,
        thickness = thickness,
        sidewallAngle = sidewallAngle_radians,
        core = Si,
        cladding=SiO2
    )

    geometry = [waveguide]

    # Set up the simulation grid
    xLocs = [0,bottomFace/2,bottomFace/2+0.1,xWidth/2]
    xVals = [dcore,dcore,dcladding,dcladding]
    xLeft = makeGrid(xLocs,xVals)
    xx    = xLeft # Only do half of the plane
    yLocs = [0,thickness/2,thickness/2+0.1,yWidth/2]
    yVals = [dcore,dcore,dcladding,dcladding]
    yLeft = makeGrid(yLocs,yVals)
    yy    = yLeft - yLeft[-1]
    yy    = np.concatenate((-(np.flip(yLeft[1:-1],0)),yLeft))

    # Run the simulation
    sim = pm.Simulation(
        geometry=geometry,
        wavelength=wavelength,
        numModes=numModes,
        xGrid=xx,
        yGrid=yy,
        background = SiO2
        )

    sim.run()
    '''
    eps = sim.getEps()
    X,Y = np.meshgrid(xx,yy)
    plt.pcolor(X,Y,np.real(eps), cmap='gray',alpha=0.5)
    plt.savefig('results.png')
    '''


    k0, Ex = sim.getFieldComponent('hr',0)
    neff = k0 / (2*np.pi/wavelength)
    return neff


# ---------------------------------------------------------------------------- #
# Run main routine
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":

    # Parse the input parameters
    kwargs         = parseInput()
    jobNumber      = kwargs['jobNumber']
    numJobs        = kwargs['numJobs']
    numModes       = kwargs['numModes']
    folderName     = kwargs['folderName']
    filenamePrefix = kwargs['filenamePrefix']
    dcore          = kwargs['dcore']
    dcladding      = kwargs['dcladding']
    xWidth         = kwargs['xWidth']
    yWidth         = kwargs['yWidth']


    # Setup the sweep
    p = {}
    p['wavelength']     = np.linspace(1.45, 1.65, 100)
    p['width']          = np.linspace(0.4, 0.6, 20)
    p['thickness']      = np.linspace(0.18, 0.24, 10)
    p['sidewallAngle']  = np.linspace(80, 90, 10)

    s = pscan.Scan.from_dict(p)

    numCombinations = s.count_total_combinations()
    blockSize = int(numCombinations / numJobs + 0.5)
    
    # Preallocate the data arrays
    neff          = np.zeros((blockSize,numModes),dtype=np.complex128)
    wavelength    = np.zeros((blockSize,))
    width         = np.zeros((blockSize,))
    thickness     = np.zeros((blockSize,))
    sidewallAngle = np.zeros((blockSize,))

    # Run the job
    start = jobNumber*blockSize
    stop  = (jobNumber + 1) * blockSize
    k     = 0
    for params in s.params(start, stop):
        # save the parameters
        wavelength[k]    = params['wavelength']
        width[k]         = params['width']
        thickness[k]     = params['thickness']
        sidewallAngle[k] = params['sidewallAngle']

        # run the iteration
        neff[k,0] = runSimulation(**params,dcore=dcore,dcladding=dcladding,xWidth=xWidth,yWidth=yWidth)
        k += 1
    print(neff)
    # Save the data
    np.savez(folderName + filenamePrefix + '_data.npz',wavelength=wavelength,width=width,thickness=thickness,sidewallAngle=sidewallAngle,neff=neff)
    print('SWEEP SUCCESSFULLY FINISHED')