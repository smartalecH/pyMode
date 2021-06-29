"""Simulation classes for front-facing access.

Alec Hammond. 2019-02-14.
"""


from subprocess import Popen, PIPE
from enum import Enum
import os
import numpy as np
import matplotlib.pyplot as plt
import pyMode as pm
import pickle

# --------------------------------------------------------------------- #
# Data structures
# --------------------------------------------------------------------- #

class Mode():
    def __init__(self,Er,Ez,Ep,Hr,Hz,Hp,k,wavelength,radius):
        self.Er = np.squeeze(Er)
        self.Ez = np.squeeze(Ez)
        self.Ep = np.squeeze(Ep)
        self.Hr = np.squeeze(Hr)
        self.Hz = np.squeeze(Hz)
        self.Hp = np.squeeze(Hp)
        self.k = k
        self.wavelength = wavelength
        self.radius=radius

        self.neff = self.wavelength * self.k / (2*np.pi)
        self.loss_um = np.imag(self.neff)*4*np.pi*4.34/self.wavelength
        self.loss_cm = self.loss_um*1e4
        self.loss_90 = self.loss_um * (2*np.pi*self.radius*0.25)
        self.Q = 2*np.pi/np.imag(self.neff)
    
    def save_mode(self,filename):
        with open(filename, 'wb') as f:
            pickle.dump(self,f)
    
    @classmethod
    def load(cls, filename):
        with open(filename, 'rb') as f:
            return pickle.load(f)
    
    def get_intensity(self):
        return self.Er**2+self.Ep**2+self.Ez**2


# --------------------------------------------------------------------- #
# Simulation Class
# --------------------------------------------------------------------- #


class Simulation:
    """Class defining parameters for a simulation with WGMS3D, as well as methods for accessing results."""
    def __init__(self, geometry, wavelength, numModes, xGrid, yGrid, radius=0, eigStart=None, boundaries=None,
                 background=pm.AIR, filenamePrefix='', folderName='',nprocs=None):
        # initialize all variables
        self.simRun = False
        self.geometry = geometry
        self.wavelength = wavelength
        self.numModes = numModes
        self.xGrid = xGrid
        self.yGrid = yGrid
        self.radius = radius
        self.eigStart = eigStart
        self.boundaries = [] if boundaries is None else boundaries
        self.background = background
        self.nprocs=nprocs

        self.filenamePrefix = filenamePrefix
        self.geomFileName = None

        # ensure folder is created (if requested)
        self.folderName = folderName
        if folderName:
            self.folderName = folderName + '/'
            if not os.path.exists(folderName):
                os.makedirs(folderName)

    def run(self):
        """Perform the simulation using WGMS3D."""
        # write the grid files
        self.makeGrid()

        # write the geometry file
        self.writeGeometry()

        # formulate the shell command
        if self.folderName == "":
            command = ""
        else:
            command = "cd " + self.folderName + " && "
        
        # add in mpi
        if self.nprocs is not None:
            command += "mpirun -np {} wgms3d".format(int(self.nprocs))
        else:
            command += "wgms3d"

        # add radius of curvature
        if self.radius:
            command += " -R {:e}".format(self.radius)

        command += " -g {}".format(self.filenamePrefix + "geometry.mgp")  # add geometry info
        command += " -l {:e}".format(self.wavelength)  # add wavelength info
        command += " -U xx.txt -V yy.txt"  # add grid info
        command += " -e -E -F -G -H"  # specify to output all fields
        command += " -n {:d}".format(self.numModes)  # specify number of modes to solve for

        if self.eigStart is not None:
            command += " -s {:e}".format(self.eigStart)  # specify initial index guess

        # Parse the boundary conditions
        for k in range(len(self.boundaries)):
            command += self.boundaries[k].output_command()

        print(command)

        # run the simulation
        process = Popen([command], stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = process.communicate()
        print("".join(chr(x) for x in stderr))
        print("".join(chr(x) for x in stdout))

        # mark the simulation as completed, allowing methods to access results
        self.simRun = True

    def getFieldComponent(self, basename, modeNumber, loadField=True):
        """Get the results for the specified field from the output folder.

        Args:
            basename (str): name of field component, e.g. 'hz'.
            modeNumber (int): number assigned to the specific mode desired.
            loadField (bool): if True, return the field component of the data. If False, only return the values
                              corresponding to wave number.
        Returns:
            (complex): wave number.
            (np.ndarray): field component requested. Returned if loadField is set.
        """
        if not self.simRun:
            raise ValueError('you must run a simulation first')

        modeNumber = "{0:0=2d}".format(modeNumber)
        filename = '{}{}-{}.bin'.format(self.folderName, basename, modeNumber)
        numX = self.xGrid.shape[0]
        numY = self.yGrid.shape[0]
        data = np.fromfile(filename, np.float64)

        # remove first two values that correspond to wavenumber
        k0 = data[0] + 1j*data[1]
        data = data[2:]

        # Parse the file
        if loadField:
            real_part_indices = np.arange(0, data.shape[0] - 1, 2)
            imag_part_indices = real_part_indices + 1
            data = data[real_part_indices] + 1j * data[imag_part_indices]
            data = data.reshape((numY, numX), order='C')

            return k0, data
        return k0

    def getFields(self):
        """Return all six fields for each mode.

        Returns:
            (list of scalars): waveNumbers
            (list of Fields): Field objects
        """
        if not self.simRun:
            raise ValueError('you must run a simulation first')

        waveNumbers = np.zeros((self.numModes,), dtype=np.complex128)

        mode_fields = []

        for mode_iter in range(self.numModes):
            modeNumber = mode_iter

            # Record data for Hr
            k0, Hr = self.getFieldComponent('hr', modeNumber)
            waveNumbers[mode_iter] = k0

            # Record data for Hz
            k0, Hz = self.getFieldComponent('hz', modeNumber)

            # Record data for Hphi
            k0, Hphi = self.getFieldComponent('hp', modeNumber)

            # Record data for Er
            k0, Er = self.getFieldComponent('er', modeNumber)

            # Record data for Ez
            k0, Ez = self.getFieldComponent('ez', modeNumber)

            # Record data for Ephi
            k0, Ephi = self.getFieldComponent('ep', modeNumber)

            mode_fields.append(Mode(Er, Ez, Ephi, Hr, Hz, Hphi, k0, self.wavelength, self.radius))
        
        return waveNumbers, mode_fields

    def getEps(self):
        """Get the eps (the geometry) for the simulation.

        Returns:
            (np.ndarray): eps (the geometry).
        """
        if not self.simRun:
            raise ValueError('you must run a simulation first')

        filename = '{}epsis.bin'.format(self.folderName)
        numX = self.xGrid.size
        numY = self.yGrid.size
        data = np.fromfile(filename, np.float64)

        # Parse the file
        real_part_indices = np.arange(0,data.shape[0] - 1, 2)
        imag_part_indices = real_part_indices + 1
        data = data[real_part_indices] + 1j * data[imag_part_indices]
        data = data.reshape((numY, numX), order='C')
        return data

    def getWavenumbers(self):
        """Get the wave numbers for the simulation.

        Returns:
            (np.ndarray): wave numbers for each mode.
        """
        if not self.simRun:
            raise ValueError('you must run a simulation first')

        waveNumbers = np.zeros((self.numModes,), dtype=np.complex128)
        for mode_iter in range(self.numModes):
            modeNumber = int("{0:0=2d}".format(mode_iter))
            k0 = self.getFieldComponent('hr', modeNumber, loadField=False)
            waveNumbers[mode_iter] = k0
        return waveNumbers

    def writeGeometry(self):
        """Initialize the geometry file for use by WGMS3D."""
        self.geomFileName = self.folderName + self.filenamePrefix + "geometry.mgp"
        print(self.background.get_n(1 / self.wavelength))
        fileContents = 'n ({:e},{:e}) \n'.format(
            np.real(self.background.get_n(1 / self.wavelength)), np.imag(self.background.get_n(1 / self.wavelength))
        )

        # Run through all the components and write them to the file
        for k in range(len(self.geometry)):
            fileContents += self.geometry[k].writeContents(self.wavelength)

        # End and write the file
        fileContents += 'x'
        with open(self.geomFileName, "w") as text_file:
            text_file.write(fileContents)

    def makeGrid(self):
        """Create the grid files for use by WGMS3D."""
        np.savetxt(self.folderName +'xx.txt', np.insert(self.xGrid, 0, int(self.xGrid.size), axis=0))
        np.savetxt(self.folderName +'yy.txt', np.insert(self.yGrid, 0, int(self.yGrid.size), axis=0))

    def plotGeometry(self, showGrid=True):
        """Plot the eps (the geometry) over the grid.
        
        Args:
            showGrid (bool): If True displays slicing grid. If False hides grids and doesn't necessarily display 'to size'. Defaults to True
            """
        eps = self.getEps()
        if showGrid:
            X,Y = np.meshgrid(self.xGrid, self.yGrid)
            plt.pcolor(X, Y, np.real(eps), cmap='gray', alpha=0.5)
        else:
            plt.imshow(np.real(eps), cmap='binary')

    def plotFields(self,modeNum=1,showGeometry=False,res=20):
        """Plots each of the field components. (Not necessarily to scale due to variability in slicing grid)

        Args:
            modeNum (int): index of mode to get. If self.numModes is 1, is basically ignored (only 1 mode has been found)
                            Must be less than self.numModes
            showGeometry (bool): If True overlays geometry on top of fields. False by default
            """

        if not self.simRun:
            raise ValueError('you must run a simulation first')
        if modeNum > self.numModes:
            raise ValueError("That mode hasn't been solved for")

        from scipy import interpolate

        fields = self.getFields()
        # Plot the fields
        titles = ['$H_r$','$H_z$','$H_{phi}$','$E_r$','$E_z$','$E_{phi}$']
        #titles = ['$H_x$','$H_y$','$H_z$','$E_x$','$E_y$','$E_z$']
        if showGeometry:
            eps = self.getEps()
        
        xnew = np.arange(self.xGrid[0],self.xGrid[-1],1/res)
        ynew = np.arange(self.yGrid[0],self.yGrid[-1],1/res)
        for k in range(6):
            #plt them, assuring that 0 is the middle value
            if self.numModes == 1:
                temp_field = np.real(np.squeeze(fields[k+1])).transpose()
            else:
                temp_field = np.real(np.squeeze(fields[k+1]))[modeNum-1,:,:].transpose()
            plt.subplot(2, 3, k+1, adjustable='box', aspect=temp_field.shape[0] / temp_field.shape[1])

            v = max( abs(temp_field.min()), abs(temp_field.max()) )

            f = interpolate.interp2d(self.xGrid, self.yGrid, np.flipud(temp_field), kind='cubic')
            
            plt.imshow(f(xnew,ynew), cmap='RdBu', vmin=-v, vmax=v)
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.axis('off')
            plt.title(titles[k])
            f = interpolate.interp2d(self.xGrid, self.yGrid, np.flipud(np.real(eps)), kind='cubic')
            if showGeometry:
                plt.imshow(f(xnew,ynew), cmap='binary', alpha=0.2)

        plt.tight_layout(pad=0.2)


# --------------------------------------------------------------------- #
# Boundary Classes
# --------------------------------------------------------------------- #


class Location(Enum):
    N = 'n'
    S = 's'
    E = 'e'
    W = 'w'


class Boundaries():
    """Base class for boundary conditions applied to simulations."""
    def __init__(self, location):
        self.location = location
    def output_command(self):
        """Return the flag syntax to specify the boundary condition to WGMS3D.

        The default is no boundary, so no flag.
        """
        command = ""
        return command


class Magnetic(Boundaries):
    """Defines magnetic boundary conditions to WGMS3D."""
    def output_command(self):
        """Return the flag syntax to specify the boundary condition to WGMS3D."""
        command = " -M {}".format(self.location.value)
        return command


class PML(Boundaries):
    """Defines PML boundary conditions to WGMS3D."""
    def __init__(self, location, thickness=2, strength=1.0):
        super(PML, self).__init__(location)
        self.thickness = thickness
        self.strength = strength

    def output_command(self):
        """Return the flag syntax to specify the boundary condition to WGMS3D."""
        command = " -P {}:{}:{}".format(self.location.value, self.thickness, self.strength)
        return command

# --------------------------------------------------------------------- #
# Utility routines
# --------------------------------------------------------------------- #

def inner_product(mode1,mode2,xgrid,ygrid):
    '''\int{E_1 \cross H^*_2 \dot dS}'''
    product = mode1.Er * np.conj(mode2.Hz) - mode1.Ez * np.conj(mode2.Hr)
    return np.trapz(np.trapz(product,xgrid,axis=1),ygrid)

def overlap_integral(mode1,mode2,xgrid,ygrid):
    prod1 = inner_product(mode1,mode2,xgrid,ygrid)
    prod2 = inner_product(mode2,mode1,xgrid,ygrid)
    norm1 = inner_product(mode1,mode1,xgrid,ygrid)
    norm2 = inner_product(mode2,mode2,xgrid,ygrid)
    return np.real(prod1*prod2/norm1)/np.real(norm2)

def fresnel_transmission(mode1,mode2):
    n1 = mode1.neff
    n2 = mode2.neff
    T = 1-np.abs((n1-n2)/(n1+n2))**2
    print("T, {}".format(T))
    return T

def power_coupling(mode1,mode2,xgrid,ygrid):
    overlap = overlap_integral(mode1,mode2,xgrid,ygrid)
    fresnel = fresnel_transmission(mode1,mode2)
    loss = overlap * fresnel
    loss_db = -10*np.log10(loss)
    return loss_db

def makeGrid(xs,hs):
    grid = [min(xs)]
    interp = interp1d(xs, hs, kind='linear')
    while grid[-1] < max(xs):
        h = interp(grid[-1])
        grid.append(grid[-1] + h)
    return np.array(grid)