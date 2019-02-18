from __future__ import division

import functools
import math
import numbers
import operator
import warnings
from collections import namedtuple
from copy import deepcopy
from numbers import Number

import numpy as np



# --------------------------------------------------------------------- #
# Geometry Classes
# --------------------------------------------------------------------- #

class Line():
    def __init__(self, nLeft, nRight, X1,Y1,X2,Y2, *args, **kwargs):
        self.nLeft = nLeft
        self.nRight = nRight
        self.X1 = X1
        self.X2 = X2
        self.Y1 = Y1
        self.Y2 = Y2
    
    def writeContents(self):
            return "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
            self.nLeft.real,self.nLeft.imag,self.nRight.real,self.nRight.imag,self.X1,self.Y1,self.X2,self.Y2
        )

class Bezier():
    def __init__(self, nLeft, nRight, X1,Y1,X2,Y2, *args, **kwargs):
        self.nLeft = nLeft
        self.nRight = nRight
        self.X1 = X1
        self.X2 = X2
        self.Y1 = Y1
        self.Y2 = Y2
    
    def writeContents(self):
            return "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
            self.nLeft.real,self.nLeft.imag,self.nRight.real,self.nRight.imag,self.X1,self.Y1,self.X2,self.Y2
        )

class Rectangle():
    def __init__(self, center, size, core, cladding, rc=0, *args, **kwargs):
        self.centerX = center.x
        self.centerY = center.y
        self.width  = size.x
        self.thickness = size.y
        self.core = core
        self.cladding = cladding
        self.rc = rc

    def writeContents(self,wavelength):
        Lines = [
            [[self.centerX-self.width/2+self.rc,self.centerY+self.thickness/2],[self.centerX+self.width/2-self.rc,self.centerY+self.thickness/2]], # top line
            [[self.centerX-self.width/2+self.rc,self.centerY-self.thickness/2],[self.centerX+self.width/2-self.rc,self.centerY-self.thickness/2]], # bottom line
            [[self.centerX-self.width/2,self.centerY-self.thickness/2+self.rc],[self.centerX-self.width/2,self.centerY+self.thickness/2-self.rc]], # left line
            [[self.centerX+self.width/2,self.centerY-self.thickness/2+self.rc],[self.centerX+self.width/2,self.centerY+self.thickness/2-self.rc]]  # right line
        ]

        Corners = [
            [[self.centerX-self.width/2,self.centerY+self.thickness/2-self.rc],[self.centerX-self.width/2,self.centerY+self.thickness/2],[self.centerX-self.width/2+self.rc,self.centerY+self.thickness/2]], # top left
            [[self.centerX-self.width/2,self.centerY-self.thickness/2+self.rc],[self.centerX-self.width/2,self.centerY-self.thickness/2],[self.centerX-self.width/2+self.rc,self.centerY-self.thickness/2]], # bottom left
            [[self.centerX+self.width/2-self.rc,self.centerY+self.thickness/2],[self.centerX+self.width/2,self.centerY+self.thickness/2],[self.centerX+self.width/2,self.centerY+self.thickness/2-self.rc]], # top right
            [[self.centerX+self.width/2-self.rc,self.centerY-self.thickness/2],[self.centerX+self.width/2,self.centerY-self.thickness/2],[self.centerX+self.width/2,self.centerY-self.thickness/2+self.rc]]  # bottom right
        ]

        cladding = self.cladding.get_n(1/wavelength)
        core     = self.core.get_n(1/wavelength)
        Indices = [
            [cladding,core],      # top line
            [core,cladding],      # bottom line
            [cladding,core],      # left line
            [core,cladding]       # right line
        ]

        fileContents = ""
        for k in range(len(Lines)):

            X1 = Lines[k][0][0]
            Y1 = Lines[k][0][1]
            X2 = Lines[k][1][0]
            Y2 = Lines[k][1][1]
            nLeft  = Indices[k][0]
            nRight = Indices[k][1]
            currentLine = "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
                nLeft.real,nLeft.imag,nRight.real,nRight.imag,X1,Y1,X2,Y2
            )
            fileContents += (currentLine)

        if self.rc > 0:
            for k in range(len(Lines)):

                X1 = Corners[k][0][0]
                Y1 = Corners[k][0][1]
                X2 = Corners[k][1][0]
                Y2 = Corners[k][1][1]
                X3 = Corners[k][2][0]
                Y3 = Corners[k][2][1]
                nLeft  = Indices[k][0]
                nRight = Indices[k][1]
                currentLine = "b ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e} {:e} {:e}\n".format(
                    nLeft.real,nLeft.imag,nRight.real,nRight.imag,X1,Y1,X2,Y2,X3,Y3
                )
                fileContents += (currentLine)
        return fileContents

class Trapezoid():
    def __init__(self, center, topFace,thickness,sidewallAngle, core, cladding, rc=0, *args, **kwargs):
        self.centerX = center.x
        self.centerY = center.y
        self.topFace  = topFace
        self.thickness = thickness
        self.sidewallAngle = sidewallAngle
        self.core = core
        self.cladding = cladding
        self.rc = rc

        self.bottomFace = topFace + 2*(thickness/np.tan(sidewallAngle))

        '''
                                      (x1,y1) ----------- (x2,y2)
                                             /            \
                                            /              \
                                           /                \
                                 (x3,y3)  --------------------  (x4,y4)
        '''

        self.X1 = self.centerX - self.topFace/2
        self.X2 = self.centerX + self.topFace/2
        self.X3 = self.centerX - self.bottomFace/2
        self.X4 = self.centerX + self.bottomFace/2
        self.Y1 = self.centerY + self.thickness / 2
        self.Y2 = self.Y1
        self.Y3 = self.centerY - self.thickness / 2
        self.Y4 = self.Y3
    def writeContents(self,wavelength):
        Lines = [
            [[self.X1 + self.rc,self.Y1],[self.X2-self.rc,self.Y2]], # top line
            [[self.X3 + self.rc,self.Y3],[self.X4-self.rc,self.Y4]], # bottom line
            [[self.X3,self.Y3 + self.rc],[self.X1,self.Y1-self.rc]], # left line
            [[self.X4,self.Y4 + self.rc],[self.X2,self.Y2-self.rc]]  # right line
        ]
        if self.rc > 0:
            print('Rounded corners not supported for trapezoids yet!')
            raise
        # TODO: Enable rounded corners
        Corners = [
            [[],[],[]], # top left
            [[],[],[]], # bottom left
            [[],[],[]], # top right
            [[],[],[]]  # bottom right
        ]

        cladding = self.cladding.get_n(1/wavelength)
        core     = self.core.get_n(1/wavelength)
        Indices = [
            [cladding,core],      # top line
            [core,cladding],      # bottom line
            [cladding,core],      # left line
            [core,cladding]       # right line
        ]

        fileContents = ""
        for k in range(len(Lines)):

            X1 = Lines[k][0][0]
            Y1 = Lines[k][0][1]
            X2 = Lines[k][1][0]
            Y2 = Lines[k][1][1]
            nLeft  = Indices[k][0]
            nRight = Indices[k][1]
            currentLine = "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
                nLeft.real,nLeft.imag,nRight.real,nRight.imag,X1,Y1,X2,Y2
            )
            fileContents += (currentLine)

        if self.rc > 0:
            for k in range(len(Lines)):

                X1 = Corners[k][0][0]
                Y1 = Corners[k][0][1]
                X2 = Corners[k][1][0]
                Y2 = Corners[k][1][1]
                X3 = Corners[k][2][0]
                Y3 = Corners[k][2][1]
                nLeft  = Indices[k][0]
                nRight = Indices[k][1]
                currentLine = "b ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e} {:e} {:e}\n".format(
                    nLeft.real,nLeft.imag,nRight.real,nRight.imag,X1,Y1,X2,Y2,X3,Y3
                )
                fileContents += (currentLine)
        return fileContents
# --------------------------------------------------------------------- #
# Material Classes
# --------------------------------------------------------------------- #

def check_nonnegative(prop, val):
    if val >= 0:
        return val
    else:
        raise ValueError("{} cannot be negative. Got {}".format(prop, val))


class Vector3(object):

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = float(x) if type(x) is int else x
        self.y = float(y) if type(y) is int else y
        self.z = float(z) if type(z) is int else z

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z

        return Vector3(x, y, z)

    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z

        return Vector3(x, y, z)

    def __mul__(self, other):
        if type(other) is Vector3:
            return self.dot(other)
        elif isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError("No operation known for 'Vector3 * {}'".format(type(other)))

    def __truediv__(self, other):
        if type(other) is Vector3:
            return Vector3(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, Number):
            return Vector3(self.x / other, self.y / other, self.z / other)
        else:
            raise TypeError("No operation known for 'Vector3 / {}'".format(type(other)))

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError("No operation known for '{} * Vector3'".format(type(other)))

    def __getitem__(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        else:
            raise IndexError("No value at index {}".format(i))

    def __repr__(self):
        return "Vector3<{}, {}, {}>".format(self.x, self.y, self.z)

    def __array__(self):
        return np.array([self.x, self.y, self.z])

    def conj(self):
        return Vector3(self.x.conjugate(), self.y.conjugate(), self.z.conjugate())

    def scale(self, s):
        x = self.x * s
        y = self.y * s
        z = self.z * s

        return Vector3(x, y, z)

    def dot(self, v):
        return self.x * v.x + self.y * v.y + self.z * v.z

    def cdot(self, v):
        return self.conj().dot(v)

    def cross(self, v):
        x = self.y * v.z - self.z * v.y
        y = self.z * v.x - self.x * v.z
        z = self.x * v.y - self.y * v.x

        return Vector3(x, y, z)

    def norm(self):
        return math.sqrt(abs(self.dot(self)))

    def unit(self):
        return self.scale(1 / self.norm())

    def close(self, v, tol=1.0e-7):
        return (abs(self.x - v.x) <= tol and
                abs(self.y - v.y) <= tol and
                abs(self.z - v.z) <= tol)

    def rotate(self, axis, theta):
        u = axis.unit()
        vpar = u.scale(u.dot(self))
        vcross = u.cross(self)
        vperp = self - vpar
        return vpar + (vperp.scale(math.cos(theta)) + vcross.scale(math.sin(theta)))

    # rotate vectors in lattice/reciprocal coords (note that the axis
    # is also given in the corresponding basis):

    def rotate_lattice(self, axis, theta, lat):
        a = lattice_to_cartesian(axis, lat)
        v = lattice_to_cartesian(self, lat)
        return cartesian_to_lattice(v.rotate(a, theta), lat)

    def rotate_reciprocal(self, axis, theta, lat):
        a = reciprocal_to_cartesian(axis, lat)
        v = reciprocal_to_cartesian(self, lat)
        return cartesian_to_reciprocal(v.rotate(a, theta), lat)


class Medium(object):

    def __init__(self, epsilon_diag=Vector3(1, 1, 1),
                 epsilon_offdiag=Vector3(0j, 0j, 0j),
                 mu_diag=Vector3(1, 1, 1),
                 mu_offdiag=Vector3(0j, 0j, 0j),
                 E_susceptibilities=[],
                 H_susceptibilities=[],
                 E_chi2_diag=Vector3(),
                 E_chi3_diag=Vector3(),
                 H_chi2_diag=Vector3(),
                 H_chi3_diag=Vector3(),
                 D_conductivity_diag=Vector3(),
                 B_conductivity_diag=Vector3(),
                 epsilon=None,
                 index=None,
                 mu=None,
                 chi2=None,
                 chi3=None,
                 D_conductivity=None,
                 B_conductivity=None,
                 E_chi2=None,
                 E_chi3=None,
                 H_chi2=None,
                 H_chi3=None,
                 valid_freq_range=None):

        if epsilon:
            epsilon_diag = Vector3(epsilon, epsilon, epsilon)
        elif index:
            i2 = index * index
            epsilon_diag = Vector3(i2, i2, i2)

        if mu:
            mu_diag = Vector3(mu, mu, mu)

        if D_conductivity:
            D_conductivity_diag = Vector3(D_conductivity, D_conductivity, D_conductivity)
        if B_conductivity:
            B_conductivity_diag = Vector3(B_conductivity, B_conductivity, B_conductivity)

        if E_chi2:
            E_chi2_diag = Vector3(E_chi2, E_chi2, E_chi2)
        if E_chi3:
            E_chi3_diag = Vector3(E_chi3, E_chi3, E_chi3)
        if H_chi2:
            H_chi2_diag = Vector3(H_chi2, H_chi2, H_chi2)
        if H_chi3:
            H_chi3_diag = Vector3(H_chi3, H_chi3, H_chi3)

        self.epsilon_diag = epsilon_diag
        self.epsilon_offdiag = epsilon_offdiag
        self.mu_diag = mu_diag
        self.mu_offdiag = mu_offdiag
        self.E_susceptibilities = E_susceptibilities
        self.H_susceptibilities = H_susceptibilities
        self.E_chi2_diag = Vector3(chi2, chi2, chi2) if chi2 else E_chi2_diag
        self.E_chi3_diag = Vector3(chi3, chi3, chi3) if chi3 else E_chi3_diag
        self.H_chi2_diag = H_chi2_diag
        self.H_chi3_diag = H_chi3_diag
        self.D_conductivity_diag = D_conductivity_diag
        self.B_conductivity_diag = B_conductivity_diag
        self.valid_freq_range = valid_freq_range

    def get_eps(self,freq):

        # Clean the input
        freq = np.squeeze([freq])

        # Compensate for scalar inputs
        if len(freq.shape) == 0:
            freq = np.array([freq])

        # Initialize with instantaneous dielectric tensor, use numpy arrays for
        # convenience and speed.
        eps = np.array([self.epsilon_diag.x,self.epsilon_diag.y,self.epsilon_diag.z],dtype='complex128')

        # Iterate through and sum up susceptibilities
        for k in range(len(self.E_susceptibilities)):
            eps = eps + self.E_susceptibilities[k].eval_susceptibility(freq)

        # Account for conductivity term
        eps = (1 + 1j*np.array(
            [self.D_conductivity_diag.x,self.D_conductivity_diag.y,self.D_conductivity_diag.z])/freq.reshape(-1,1)) * eps
        eps = [Vector3(eps[row,0],eps[row,1],eps[row,2]) for row in range(freq.shape[0])]
        
        # Only return x-component for now
        return np.squeeze(np.array([o.x for o in eps]))
    
    def get_n(self,freq):
        return np.sqrt(self.get_eps(freq))


class Susceptibility(object):

    def __init__(self, sigma_diag=Vector3(), sigma_offdiag=Vector3(), sigma=None):
        self.sigma_diag = Vector3(sigma, sigma, sigma) if sigma else sigma_diag
        self.sigma_offdiag = sigma_offdiag

    def transform(self, m):
        sigma = Matrix(Vector3(self.sigma_diag.x, self.sigma_offdiag.x, self.sigma_offdiag.y),
                       Vector3(self.sigma_offdiag.x, self.sigma_diag.y, self.sigma_offdiag.z),
                       Vector3(self.sigma_offdiag.y, self.sigma_offdiag.z, self.sigma_diag.z))
        new_sigma = (m * sigma * m.transpose()) / abs(m.determinant())
        self.sigma_diag = Vector3(new_sigma.c1.x, new_sigma.c2.y, new_sigma.c3.z)
        self.sigma_offdiag = Vector3(new_sigma.c2.x, new_sigma.c3.x, new_sigma.c3.y)


class LorentzianSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super(LorentzianSusceptibility, self).__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma

    def eval_susceptibility(self,freq):

        # Use numpy arrays for element-wise multiplication later
        sigma = np.array([self.sigma_diag.x,self.sigma_diag.y,self.sigma_diag.z])
        eps = np.zeros((freq.shape[0],3),dtype='complex128')

        # Iterate through dimensions
        for k in range(3):
            eps[:,k] = (sigma[k]*self.frequency*self.frequency) / (self.frequency*self.frequency - freq*freq - 1j*self.gamma*freq)
        return eps


class DrudeSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super(DrudeSusceptibility, self).__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma

    def eval_susceptibility(self,freq):

        # Use numpy arrays for element-wise multiplication later
        sigma = np.array([self.sigma_diag.x,self.sigma_diag.y,self.sigma_diag.z])
        eps = np.zeros((freq.shape[0],3),dtype='complex128')

        # Iterate through dimensions
        for k in range(3):
            eps[:,k] = (-sigma[k]*self.frequency*self.frequency) / (freq*(freq + 1j*self.gamma))
        return eps