import subprocess
from subprocess import Popen, PIPE
import json
import numpy as np
import os
from enum import Enum
import pyMode as pm

# --------------------------------------------------------------------- #
# Simulation Class
# --------------------------------------------------------------------- #

class Simulation:
    def __init__(self, geometry, wavelength, numModes, xGrid, yGrid, radius, eigStart=None, boundaries=[], background = pm.AIR, filenamePrefix = '', folderName = '', *args, **kwargs):
        # initialize all variables
        self.simRun = False
        self.geometry = geometry
        self.wavelength = wavelength
        self.numModes = numModes
        self.xGrid = xGrid
        self.yGrid = yGrid
        self.radius = radius
        self.eigStart = eigStart
        self.boundaries = boundaries
        self.background = background
        

        self.filenamePrefix = filenamePrefix

        # ensure folder is created (if requested)
        self.folderName = folderName
        if folderName is not "":
            self.folderName = folderName + '/'
            if not os.path.exists(folderName):
                os.makedirs(folderName)
        

        return
    def run(self):

        # write the grid files
        self.makeGrid()

        # write the geometry file
        self.writeGeometry()

        # formulate the shell command
        if self.folderName is "":
            command = "wgms3d"
        else:
            command  = "cd " + self.folderName + " && wgms3d"    

        # add radius of curvature
        if self.radius is not 0:
            command += "R {:e}.format(self.radius)"

        command += " -g {}".format(self.geomFileName)            # add geometry info
        command += " -l {:e}".format(self.wavelength)            # add wavelength info
        command += " -U xx.txt -V yy.txt"                        # add grid info
        command += " -e -E -F -G -H"                             # specify to output all fields
        command += " -n {:d}".format(self.numModes)              # specify number of modes to solve for
        
        if self.eigStart is not None:
            command += " -s {:e}".format(self.eigStart)          # specify initial index guess
        
        # Parse the boundary conditions
        for k in range(len(self.boundaries)):
            command += self.boundaries[k].output_command()

        print(command)

        # run the simulation
        process = Popen([command], stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = process.communicate()
        print("".join( chr(x) for x in stderr))
        print("".join( chr(x) for x in stdout))

        self.simRun = True
        return
    def getFieldComponent(self,basename,modeNumber):
        if self.simRun:

            filename = '{}{}-{}.bin'.format(self.folderName, basename, modeNumber)
            numX = self.xGrid.shape[0]
            numY = self.yGrid.shape[0]
            data = np.fromfile(filename, np.float64)

            # remove first two values that correspond to wavenumber
            k0 = data[0] + 1j*data[1]
            data = data[2:]

            # Parse the file
            real_part_indices = np.arange(0,data.shape[0]-1,2)
            imag_part_indices = real_part_indices + 1
            data = data[real_part_indices] + 1j*data[imag_part_indices]
            data = (data.reshape((numY,numX)))

            return k0, data
        else:
            print('You must run a simulation first!')
            raise
    def getFields(self):
        if self.simRun:
            waveNumbers = np.zeros((self.numModes,),dtype=np.complex128)

            Hr          = np.zeros((self.numModes,self.xGrid.size,self.yGrid.size),dtype=np.complex128)
            Hz          = np.zeros((self.numModes,self.xGrid.size,self.yGrid.size),dtype=np.complex128)
            Hphi        = np.zeros((self.numModes,self.xGrid.size,self.yGrid.size),dtype=np.complex128)

            Er          = np.zeros((self.numModes,self.xGrid.size,self.yGrid.size),dtype=np.complex128)
            Ez          = np.zeros((self.numModes,self.xGrid.size,self.yGrid.size),dtype=np.complex128)
            Ephi        = np.zeros((self.numModes,self.xGrid.size,self.yGrid.size),dtype=np.complex128)

            Eps        = self.getEps().T
            k0_modes = np.zeros((self.numModes,),dtype=np.complex128)

            for mode_iter in range(self.numModes):
                modeNumber = "{0:0=2d}".format(mode_iter)
                # Record data for Hr
                k0, data = self.getFieldComponent('hr', modeNumber)
                waveNumbers[mode_iter] = k0
                Hr[mode_iter,:,:] = data.T

                # Record data for Hz
                k0, data = self.getFieldComponent('hz', modeNumber)
                Hz[mode_iter,:,:] = data.T

                # Record data for Hphi
                k0, data = self.getFieldComponent('hp', modeNumber)
                Hphi[mode_iter,:,:] = data.T

                # Record data for Er
                k0, data = self.getFieldComponent('er', modeNumber)
                Er[mode_iter,:,:] = data.T

                # Record data for Ez
                k0, data = self.getFieldComponent('ez', modeNumber)
                Ez[mode_iter,:,:] = data.T

                # Record data for Ephi
                k0, data = self.getFieldComponent('ep', modeNumber)
                Ephi[mode_iter,:,:] = data.T
                k0_modes[mode_iter] = k0

                return waveNumbers, Hr, Hz, Hphi, Er, Ez, Ephi
        else:
            print('You must run a simulation first!')
            raise
    def getEps(self):
        if self.simRun:
            filename = '{}epsis.bin'.format(directory)
            numX = self.xGrid.shape[0]
            numY = self.yGrid.shape[0]
            data = np.fromfile(filename, np.float64)

            # Parse the file
            real_part_indices = np.arange(0,data.shape[0]-1,2)
            imag_part_indices = real_part_indices + 1
            data = data[real_part_indices] + 1j*data[imag_part_indices]
            data = (data.reshape((numY,numX)))
            return data
        else:
            print('You must run a simulation first!')
            raise
    def getWavenumbers(self):
        if self.simRun:
            waveNumbers = np.zeros((self.numModes,),dtype=np.complex128)
            for mode_iter in range(self.numModes):
                modeNumber = "{0:0=2d}".format(mode_iter)
                k0, data = self.getFieldComponent('hr', modeNumber)
                waveNumbers[mode_iter] = k0
            return waveNumbers
        else:
            print('You must run a simulation first!')
            raise
    def writeGeometry(self):

        # Initialize the geometry file
        self.geomFileName = self.folderName + self.filenamePrefix + "geometry.mgp"
        print(self.background.get_n(1/self.wavelength))
        fileContents = 'n ({:e},{:e}) \n'.format(np.real(self.background.get_n(1/self.wavelength)),np.imag(self.background.get_n(1/self.wavelength)))

        # Run through all the components and write them to the file
        for k in range(len(self.geometry)):
            fileContents += self.geometry[k].writeContents(self.wavelength)
        
        # End and write the file
        fileContents += 'x'
        with open(self.geomFileName, "w") as text_file:
            text_file.write(fileContents)
    def makeGrid(self):
        np.savetxt(self.folderName +'xx.txt',np.insert(self.xGrid, 0, int(self.xGrid.size), axis=0))
        np.savetxt(self.folderName +'yy.txt',np.insert(self.yGrid, 0, int(self.yGrid.size), axis=0))
    
# --------------------------------------------------------------------- #
# Boundary Classes
# --------------------------------------------------------------------- #
class Locations(Enum):
    N = 'n'
    S = 's'
    E = 'e'
    W = 'w'

class Boundaries():
    def __init__(self, location, *args, **kwargs):
        self.location = location
    def output_command(self):
        command = ""
        return command

class Magnetic(Boundaries):
    def __init__(self, location, *args, **kwargs):
        self.location = location
    def output_command(self):
        command = " -M {}".format(self.location.value)
        return command

class PML(Boundaries):
    def __init__(self, location, thickness=2, strength=1.0, *args, **kwargs):
        self.location = location
        self.thickness = thickness
        self.strength = strength
    def output_command(self):
        command = " -P {}:{}:{}".format(self.location.value,self.thickness,self.strength)
        return command