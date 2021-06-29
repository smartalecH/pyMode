"""Definitions for waveguide geometries to be used in simulation.

Alec Hammond. 2019-02-16.
"""

from __future__ import division

import math
from numbers import Number

import numpy as np


# --------------------------------------------------------------------- #
# Geometry Classes
# --------------------------------------------------------------------- #


class Line():
    """A line defined in 2D space by start and end points."""
    def __init__(self, nLeft, nRight, X1, Y1, X2, Y2):
        """Initialize the line object.

        Args:
            nLeft (???): ???
            nRight (???): ???
            X1 (float): position of first endpoint in the x direction.
            Y1 (float): position of first endpoint in the y direction.
            X2 (float): position of second endpoint in the x direction.
            Y2 (float): position of second endpoint in the y direction.
        """
        self.nLeft = nLeft
        self.nRight = nRight
        self.X1 = X1
        self.X2 = X2
        self.Y1 = Y1
        self.Y2 = Y2

    def writeContents(self):
        """Return a string defining the line in the way understood by WGMS3D.

        Returns:
            (str): string representation of line object.
        """
        return "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
            self.nLeft.real, self.nLeft.imag, self.nRight.real, self.nRight.imag, self.X1, self.Y1, self.X2, self.Y2
        )


class Bezier():
    """A Bezier curve (???) defined in 2D space by start and end points."""
    def __init__(self, nLeft, nRight, X1, Y1, X2, Y2):
        """Initialize the curve object.

        Args:
            nLeft (???): ???
            nRight (???): ???
            X1 (Number): position of first endpoint in the x direction.
            Y1 (Number): position of first endpoint in the y direction.
            X2 (Number): position of second endpoint in the x direction.
            Y2 (Number): position of second endpoint in the y direction.
        """
        self.nLeft = nLeft
        self.nRight = nRight
        self.X1 = X1
        self.X2 = X2
        self.Y1 = Y1
        self.Y2 = Y2

    def writeContents(self):
        """Return a string defining the curve in the way understood by WGMS3D.

        Returns:
            (str): string representation of Bezier curve object.
        """
        return "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
            self.nLeft.real, self.nLeft.imag, self.nRight.real, self.nRight.imag, self.X1, self.Y1, self.X2, self.Y2
        )

class Rectangle():
    """Rectangle defined in 2D space, with the material of the core and cladding specified."""
    def __init__(self, center, size, core, cladding, rc=0):
        """Initialize the curve object.

        Args:
            center (coordinate): location of center of rectangle in 2D space; must contain x and y attributes.
            size (coordinate): width (x direction) and thickness (y direction) of rectangle; must contain x and y
                               attributes.
            core (pyMode.materials.*): material of core.
            cladding (pyMode.materials.*): material of cladding.
            rc (???): roundness of corners.
        """
        self.centerX = center.x
        self.centerY = center.y
        self.width = size.x
        self.thickness = size.y
        self.core = core
        self.cladding = cladding
        self.rc = rc

    def writeContents(self, wavelength):
        """Return a string defining the rectangle in the way understood by WGMS3D.

        Returns:
            (str): string representation of rectangle object.
        """
        # move attributes to local namespace for readability
        centerX = self.centerX
        centerY = self.centerY
        width = self.width
        thickness = self.thickness
        core = self.core
        cladding = self.cladding
        rc = self.rc

        lines = [  # top, bottom, left, right
            [[centerX - width / 2 + rc, centerY + thickness / 2], [centerX + width / 2 - rc, centerY + thickness / 2]],
            [[centerX - width / 2 + rc, centerY - thickness / 2], [centerX + width / 2 - rc, centerY - thickness / 2]],
            [[centerX - width / 2, centerY - thickness / 2 + rc], [centerX - width / 2, centerY + thickness / 2 - rc]],
            [[centerX + width / 2, centerY - thickness / 2 + rc], [centerX + width / 2, centerY + thickness / 2 - rc]]
        ]

        corners = [
            [# top-left
                [centerX - width / 2, centerY + thickness / 2 - rc],
                [centerX - width / 2, centerY + thickness / 2],
                [centerX - width / 2 + rc, centerY + thickness / 2]
            ],
            [  # bottom-left
                [centerX - width / 2, centerY - thickness / 2 + rc],
                [centerX - width / 2, centerY - thickness / 2],
                [centerX - width / 2 + rc, centerY - thickness / 2]
            ],
            [  # top-right
                [centerX + width / 2 - rc, centerY + thickness / 2],
                [centerX + width / 2, centerY + thickness / 2],
                [centerX + width / 2, centerY + thickness / 2 - rc]
            ],
            [  # bottom-right
                [centerX + width / 2 - rc, centerY - thickness / 2],
                [centerX + width / 2, centerY - thickness / 2],
                [centerX + width / 2, centerY - thickness / 2 + rc]
            ]
        ]

        cladding = cladding.get_n(1 / wavelength)
        core = core.get_n(1 / wavelength)
        indices = [
            [cladding, core],  # top line
            [core, cladding],  # bottom line
            [cladding, core],  # left line
            [core, cladding]  # right line
        ]

        fileContents = ""

        for k, line in enumerate(lines):
            X1 = line[0][0]
            Y1 = line[0][1]
            X2 = line[1][0]
            Y2 = line[1][1]
            nLeft = indices[k][0]
            nRight = indices[k][1]
            currentLine = "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
                nLeft.real, nLeft.imag, nRight.real, nRight.imag, X1, Y1, X2, Y2
            )
            fileContents += (currentLine)

        if rc > 0:
            for k in range(len(lines)):
                X1 = corners[k][0][0]
                Y1 = corners[k][0][1]
                X2 = corners[k][1][0]
                Y2 = corners[k][1][1]
                X3 = corners[k][2][0]
                Y3 = corners[k][2][1]
                nLeft = indices[k][0]
                nRight = indices[k][1]
                currentLine = "b ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e} {:e} {:e}\n".format(
                    nLeft.real, nLeft.imag, nRight.real, nRight.imag, X1, Y1, X2, Y2, X3, Y3
                )
                fileContents += (currentLine)
        return fileContents


class Trapezoid():
    r"""Trapezoid defined in 2D space, with the material of the core and cladding specified.

                                      (x1,y1) ----------- (x2,y2)
                                             /            \
                                            /              \
                                           /                \
                                 (x3,y3)  --------------------  (x4,y4)
    """
    def __init__(self, center, top_face, thickness, sidewall_angle, core, cladding, rc=0):
        """Initialize the trapezoid object.

        Args:
            center (coordinate): location of center of trapezoid in 2D space; must contain x and y attributes.
            top_face (Number): width of top face of rectangle.
            thickness (Number): thickness (y direction) of trapezoid.
            sidewall_angle (Number): angle of sidewall in radians. 0.5 * pi makes a rectangle.
            core (pyMode.materials.*): material of core.
            cladding (pyMode.materials.*): material of cladding.
            rc (???): roundness of corners.
        """
        self.centerX = center.x
        self.centerY = center.y
        self.top_face = top_face
        self.thickness = thickness
        self.sidewallAngle = sidewall_angle
        self.core = core
        self.cladding = cladding
        self.rc = rc

        # calculate length of bottom face based on the sidewall angle
        self.bottom_face = top_face + 2 * (thickness / np.tan(sidewall_angle))

        self.X1 = self.centerX - self.top_face/2
        self.X2 = self.centerX + self.top_face/2
        self.X3 = self.centerX - self.bottom_face/2
        self.X4 = self.centerX + self.bottom_face/2
        self.Y1 = self.centerY + self.thickness / 2
        self.Y2 = self.Y1
        self.Y3 = self.centerY - self.thickness / 2
        self.Y4 = self.Y3

    def writeContents(self, wavelength):
        """Return a string defining the trapezoid in the way understood by WGMS3D.

        Returns:
            (str): string representation of trapezoid object.
        """
        lines = [
            [[self.X1 + self.rc, self.Y1], [self.X2 - self.rc, self.Y2]],  # top line
            [[self.X3 + self.rc, self.Y3], [self.X4 - self.rc, self.Y4]],  # bottom line
            [[self.X3, self.Y3 + self.rc], [self.X1, self.Y1 - self.rc]],  # left line
            [[self.X4, self.Y4 + self.rc], [self.X2, self.Y2 - self.rc]]  # right line
        ]

        if self.rc > 0:
            raise ValueError('rounded corners not supported for trapezoids yet')

        # corners = [
        #     [[], [], []],  # top left
        #     [[], [], []],  # bottom left
        #     [[], [], []],  # top right
        #     [[], [], []]  # bottom right
        # ]

        cladding = self.cladding.get_n(1 / wavelength)
        core = self.core.get_n(1 / wavelength)
        indices = [
            [cladding, core],  # top line
            [core, cladding],  # bottom line
            [cladding, core],  # left line
            [core, cladding]  # right line
        ]

        fileContents = ""

        for k, line in enumerate(lines):
            X1 = line[0][0]
            Y1 = line[0][1]
            X2 = line[1][0]
            Y2 = line[1][1]
            nLeft = indices[k][0]
            nRight = indices[k][1]
            currentLine = "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
                nLeft.real, nLeft.imag, nRight.real, nRight.imag, X1, Y1, X2, Y2
            )
            fileContents += (currentLine)

        # TODO: Enable rounded corners

        return fileContents

class RibWaveguide():
    r"""Trapezoid defined in 2D space, with the material of the core and cladding specified.

                                                 width
                                      (x1,y1) ----------- (x2,y2)      |                         |
                                             /            \            | etch_depth              |
                                            /              \  (angle)  |                         |
    (x7,y7) _______________________________/                \________________________  (x4,y4)   | thickness
            |                          (x8,y8)            (x3,y3)                   |            |
            |                                                                       |            |
            |                                                                       |            |
    (x6,y6) -------------------------------------------------------------------------  (x5,y5)   |
                                                base

    """
    def __init__(self, center, width, thickness, sidewall_angle, etch_depth, base, core, cladding, rc=0):
        """Initialize the trapezoid object.

        Args:
            center (coordinate): location of center of trapezoid in 2D space; must contain x and y attributes.
            width (Number): width of top face of the trapezoid.
            thickness (Number): thickness (y direction) of entire struct.
            sidewall_angle (Number): angle of sidewall in radians. 0.5 * pi makes a rectangle.
            etch_depth (Number): height of trapezoid
            base (Number): width of slab section 
            core (pyMode.materials.*): material of core.
            cladding (pyMode.materials.*): material of cladding.
            rc (???): roundness of corners.
        """
        self.centerX = center.x
        self.centerY = center.y
        self.width = width
        self.thickness = thickness
        self.sidewallAngle = sidewall_angle
        self.etch_depth = etch_depth
        self.base = base
        self.core = core
        self.cladding = cladding
        self.rc = rc
        self.sidewall_angle=sidewall_angle

        self._compute_vertices()
    
    def _compute_vertices(self):

        # calculate length of bottom face based on the sidewall angle
        self.bottom_face = self.width + 2 * (self.etch_depth / np.tan(self.sidewall_angle))

        self.X1 = self.centerX - self.width/2
        self.X2 = self.centerX + self.width/2
        self.X3 = self.centerX + self.bottom_face/2
        self.X4 = self.centerX + self.base/2
        self.X5 = self.X4
        self.X6 = self.centerX - self.base/2
        self.X7 = self.X6
        self.X8 = self.centerX - self.bottom_face/2
        
        self.Y1 = self.centerY + self.thickness / 2
        self.Y2 = self.Y1
        self.Y3 = self.Y1 - self.etch_depth
        self.Y4 = self.Y3
        self.Y5 = self.centerY - self.thickness / 2
        self.Y6 = self.Y5
        self.Y7 = self.Y3
        self.Y8 = self.Y7

    def writeContents(self, wavelength):
        """Return a string defining the trapezoid in the way understood by WGMS3D.

        Returns:
            (str): string representation of trapezoid object.
        """
        lines = [
            [[self.X1, self.Y1], [self.X2, self.Y2]],
            [[self.X2, self.Y2], [self.X3, self.Y3]],
            [[self.X3, self.Y3], [self.X4, self.Y4]],
            #[[self.X4, self.Y4], [self.X5, self.Y5]],
            [[self.X5, self.Y5], [self.X6, self.Y6]],
            #[[self.X6, self.Y6], [self.X7, self.Y7]],
            [[self.X7, self.Y7], [self.X8, self.Y8]],
            [[self.X8, self.Y8], [self.X1, self.Y1]]
        ]

        if self.rc > 0:
            raise ValueError('rounded corners not supported for rib waveguides yet')

        # corners = [
        #     [[], [], []],  # top left
        #     [[], [], []],  # bottom left
        #     [[], [], []],  # top right
        #     [[], [], []]  # bottom right
        # ]

        cladding = self.cladding.get_n(1 / wavelength)
        core = self.core.get_n(1 / wavelength)
        indices = [
            [cladding, core]
        ]*8

        fileContents = ""

        for k, line in enumerate(lines):
            X1 = line[0][0]
            Y1 = line[0][1]
            X2 = line[1][0]
            Y2 = line[1][1]
            nLeft = indices[k][0]
            nRight = indices[k][1]
            currentLine = "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
                nLeft.real, nLeft.imag, nRight.real, nRight.imag, X1, Y1, X2, Y2
            )
            fileContents += (currentLine)

        # TODO: Enable rounded corners

        return fileContents
# --------------------------------------------------------------------- #
# Material Classes
# --------------------------------------------------------------------- #


def check_nonnegative(prop, val):
    """Ensure that val is greater than or equal to zero. Otherwise, raise a ValueError.

    Args:
        prop: property whose value must not be negative.
        val (int/float): value of prop.
    Returns:
        (int/float): val if greater than or equal to zero.
    Raises:
        ValueError: if val < 0.
    """
    if val < 0:
        raise ValueError("{} cannot be negative. Got {}".format(prop, val))
    return val


class Vector3(object):
    """Vector in 3D space."""
    def __init__(self, x=0.0, y=0.0, z=0.0):
        """Initialize the vector with coordinate (x, y, z).

        Args:
            x (Number): magnitude in x direction.
            y (Number): magnitude in y direction.
            z (Number): magnitude in z direction.
        """
        self.x = float(x) if isinstance(x, int) else x
        self.y = float(y) if isinstance(y, int) else y
        self.z = float(z) if isinstance(z, int) else z

    def __eq__(self, other):
        """Check equality.

        Args:
            other (Vector3): other vector to compare to.
        Returns:
            (bool): equality of self and other.
        """
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __add__(self, other):
        """Return the vector sum of two vectors.

        Args:
            other (Vector3): other vector to add.
        Returns:
            (Vector3): vector sum.
        """
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z

        return Vector3(x, y, z)

    def __sub__(self, other):
        """Return the vector subtraction of two vectors.

        Args:
            other (Vector3): the vector to subtract.
        Returns:
            (Vector3): vector result.
        """
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z

        return Vector3(x, y, z)

    def __mul__(self, other):
        """Infer whether vector dot product or scalar product is desired, based on type.

        Args:
            other (Vector3/Number): multiplicand.
        Returns:
            (Vector3): suitably multiplied vector result.
        """
        if isinstance(other, Vector3):
            return self.dot(other)
        elif isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError("no operation known for 'Vector3 * {}'".format(type(other)))

    def __truediv__(self, other):
        """Suitably broadcast the division based on type.

        Args:
            other (Vector3/Number): divisor.
        Returns:
            (Vector3): suitably divided vector result.
        """
        if isinstance(other, Vector3):
            return Vector3(self.x / other.x, self.y / other.y, self.z / other.z)
        elif isinstance(other, Number):
            return Vector3(self.x / other, self.y / other, self.z / other)
        else:
            raise TypeError("no operation known for 'Vector3 / {}'".format(type(other)))

    def __rmul__(self, other):
        """Broadcast the multiplication by a number to each element of the vector.

        Args:
            other (Number): multiplicand.
        Returns:
            (Vector3): suitably multiplied vector result.
        """
        if isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError("no operation known for '{} * Vector3'".format(type(other)))

    def __getitem__(self, i):
        """Order the dimensions and allow indexing for values.

        Args:
            i (int): index.
        Returns:
            (float): value at specified index.
        Raises:
            IndexError: if an index not in {0, 1, 2} is specified.
        """
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        raise IndexError("no value for index {}".format(i))

    def __repr__(self):
        """Create a string representation of the vector.

        Returns:
            (str): string representation."""
        return "Vector3<{}, {}, {}>".format(self.x, self.y, self.z)

    def __array__(self):
        """Create a NumPy array from the vector.

        Returns:
            (np.ndarray): NumPy array of vector values.
        """
        return np.array([self.x, self.y, self.z])

    def conj(self):
        """Return the complex conjugate of the vector, elementwise.

        Returns:
            (Vector3): complex conjugate of the vector.
        """
        return Vector3(self.x.conjugate(), self.y.conjugate(), self.z.conjugate())

    def scale(self, s):
        """Scale each dimension by a constant value (i.e. scalar multiplication).

        Args:
            s (Number): scalar to multiply.
        Returns:
            (Vector3): scaled vector.
        """
        x = self.x * s
        y = self.y * s
        z = self.z * s

        return Vector3(x, y, z)

    def dot(self, other):
        """Compute the dot product of the vector with another vector.

        Args:
            other (Vector3): other vector to dot product.
        Returns:
            (float): value of dot product.
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cdot(self, other):
        """Compute the conjugate dot product of the vector with another.

        Args:
            other (Vector3): other vector to dot product.
        Returns:
            (float): value of conjugate dot product.
        """
        return self.conj().dot(other)

    def cross(self, other):
        """Cross product of two vectors.

        Args:
            other (Vector3): other vector for cross product.
        Returns:
            (Vector3): cross product result.
        """
        x = self.y * other.z - self.z * other.y
        y = self.z * other.x - self.x * other.z
        z = self.x * other.y - self.y * other.x

        return Vector3(x, y, z)

    def norm(self):
        """Calculate the norm of the vector.

        Returns:
            (float): norm of vector.
        """
        return math.sqrt(abs(self.dot(self)))

    def unit(self):
        """Find the unit vector pointed in the same direction as this vector.

        Returns:
            (Vector3): normalized vector.
        """
        return self.scale(1 / self.norm())

    def close(self, other, tol=1.0e-7):
        """Determine whether two vectors are close to each other within a tolerance along all three dimensions."""
        return (abs(self.x - other.x) <= tol and
                abs(self.y - other.y) <= tol and
                abs(self.z - other.z) <= tol)

    def rotate(self, axis, theta):
        """Rotate the vector across an axis by angle theta.

        Args:
            axis (Vector3): axis of rotation.
            theta (float): angle of rotation (in radians).
        Returns:
            (Vector3): rotated vector.
        """
        u = axis.unit()
        v_par = u.scale(u.dot(self))
        v_cross = u.cross(self)
        v_perp = self - v_par
        return v_par + (v_perp.scale(math.cos(theta)) + v_cross.scale(math.sin(theta)))

    # rotate vectors in lattice/reciprocal coords (note that the axis
    # is also given in the corresponding basis):

    # def rotate_lattice(self, axis, theta, lat):
    #     a = lattice_to_cartesian(axis, lat)
    #     v = lattice_to_cartesian(self, lat)
    #     return cartesian_to_lattice(v.rotate(a, theta), lat)

    # def rotate_reciprocal(self, axis, theta, lat):
    #     a = reciprocal_to_cartesian(axis, lat)
    #     v = reciprocal_to_cartesian(self, lat)
    #     return cartesian_to_reciprocal(v.rotate(a, theta), lat)


class Medium(object):
    def __init__(self,
                 epsilon_diag=Vector3(1, 1, 1),
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

    def get_eps(self, freq):

        # Clean the input
        freq = np.squeeze([freq])

        # Compensate for scalar inputs
        if not freq.shape:
            freq = np.array([freq])

        # Initialize with instantaneous dielectric tensor, use numpy arrays for convenience and speed
        eps = np.array([self.epsilon_diag.x, self.epsilon_diag.y, self.epsilon_diag.z], dtype='complex128')

        # Iterate through and sum up susceptibilities
        for k in range(len(self.E_susceptibilities)):
            eps = eps + self.E_susceptibilities[k].eval_susceptibility(freq)

        # Account for conductivity term
        eps = (1 + 1j * np.array(
            [self.D_conductivity_diag.x, self.D_conductivity_diag.y, self.D_conductivity_diag.z]
        ) / freq.reshape(-1, 1)) * eps
        eps = [Vector3(eps[row, 0], eps[row, 1], eps[row, 2]) for row in range(freq.shape[0])]

        # Only return x-component for now
        return np.squeeze(np.array([o.x for o in eps]))

    def get_n(self, freq):
        return np.sqrt(self.get_eps(freq))


class Susceptibility(object):

    def __init__(self, sigma_diag=Vector3(), sigma_offdiag=Vector3(), sigma=None):
        self.sigma_diag = Vector3(sigma, sigma, sigma) if sigma else sigma_diag
        self.sigma_offdiag = sigma_offdiag

    def transform(self, m):
        sigma = Matrix(
            Vector3(self.sigma_diag.x, self.sigma_offdiag.x, self.sigma_offdiag.y),
            Vector3(self.sigma_offdiag.x, self.sigma_diag.y, self.sigma_offdiag.z),
            Vector3(self.sigma_offdiag.y, self.sigma_offdiag.z, self.sigma_diag.z)
                      )
        new_sigma = (m * sigma * m.transpose()) / abs(m.determinant())
        self.sigma_diag = Vector3(new_sigma.c1.x, new_sigma.c2.y, new_sigma.c3.z)
        self.sigma_offdiag = Vector3(new_sigma.c2.x, new_sigma.c3.x, new_sigma.c3.y)


class LorentzianSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super(LorentzianSusceptibility, self).__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma

    def eval_susceptibility(self, freq):

        # Use numpy arrays for element-wise multiplication later
        sigma = np.array([self.sigma_diag.x, self.sigma_diag.y, self.sigma_diag.z])
        eps = np.zeros((freq.shape[0], 3), dtype='complex128')

        # Iterate through dimensions
        for k in range(3):
            num = sigma[k] * self.frequency * self.frequency
            denom = self.frequency * self.frequency - freq * freq - 1j * self.gamma * freq
            eps[:, k] = num / denom
        return eps


class DrudeSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super(DrudeSusceptibility, self).__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma

    def eval_susceptibility(self, freq):
        # Use numpy arrays for element-wise multiplication later
        sigma = np.array([self.sigma_diag.x, self.sigma_diag.y, self.sigma_diag.z])
        eps = np.zeros((freq.shape[0], 3), dtype='complex128')

        # Iterate through dimensions
        for k in range(3):
            eps[:, k] = (-sigma[k] * self.frequency * self.frequency) / (freq * (freq + 1j * self.gamma))
        return eps
