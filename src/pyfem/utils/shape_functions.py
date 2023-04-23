from math import sqrt

from numpy import dot, empty, zeros, cross
from scipy.linalg import norm, det, inv
from scipy.special.orthogonal import p_roots as gauss_scheme


class ShapeData:
    def __init__(self):
        self.h = None
        self.dhdxi = None
        self.xi = None


class ElementShapeData:

    def __init__(self):
        self.shape_data = []

    def __iter__(self):
        return iter(self.shape_data)

    def __len__(self):
        return len(self.shape_data)


def get_shape_line2(xi):
    # Check the dimensions of the physical space
    if type(xi) != float:
        raise NotImplementedError('1D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(2)
    shape_data.dhdxi = empty(shape=(2, 1))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = 0.5 * (1.0 - xi)
    shape_data.h[1] = 0.5 * (1.0 + xi)

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -0.5
    shape_data.dhdxi[1, 0] = 0.5

    return shape_data


def get_shape_line3(xi):
    # Check the dimension of physical space
    if type(xi) != float:
        raise NotImplementedError('1D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(3)
    shape_data.dhdxi = empty(shape=(1, 3))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = 0.5 * (1.0 - xi) - 0.5 * (1.0 - xi * xi)
    shape_data.h[1] = 1 - xi * xi
    shape_data.h[2] = 0.5 * (1.0 + xi) - 0.5 * (1.0 - xi * xi)

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -0.5 + xi
    shape_data.dhdxi[0, 1] = -2.0 * xi
    shape_data.dhdxi[0, 2] = 0.5 + xi

    return shape_data


def get_shape_tria3(xi):
    # Check the dimension of physical space
    if len(xi) != 2:
        raise NotImplementedError('2D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(3)
    shape_data.dhdxi = empty(shape=(3, 2))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = 1.0 - xi[0] - xi[1]
    shape_data.h[1] = xi[0]
    shape_data.h[2] = xi[1]

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -1.0
    shape_data.dhdxi[1, 0] = 1.0
    shape_data.dhdxi[2, 0] = 0.0

    shape_data.dhdxi[0, 1] = -1.0
    shape_data.dhdxi[1, 1] = 0.0
    shape_data.dhdxi[2, 1] = 1.0

    return shape_data


def get_shape_quad4(xi):
    # Check the dimension of physical space
    if len(xi) != 2:
        raise NotImplementedError('2D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(4)
    shape_data.dhdxi = empty(shape=(4, 2))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = 0.25 * (1.0 - xi[0]) * (1.0 - xi[1])
    shape_data.h[1] = 0.25 * (1.0 + xi[0]) * (1.0 - xi[1])
    shape_data.h[2] = 0.25 * (1.0 + xi[0]) * (1.0 + xi[1])
    shape_data.h[3] = 0.25 * (1.0 - xi[0]) * (1.0 + xi[1])

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -0.25 * (1.0 - xi[1])
    shape_data.dhdxi[1, 0] = 0.25 * (1.0 - xi[1])
    shape_data.dhdxi[2, 0] = 0.25 * (1.0 + xi[1])
    shape_data.dhdxi[3, 0] = -0.25 * (1.0 + xi[1])

    shape_data.dhdxi[0, 1] = -0.25 * (1.0 - xi[0])
    shape_data.dhdxi[1, 1] = -0.25 * (1.0 + xi[0])
    shape_data.dhdxi[2, 1] = 0.25 * (1.0 + xi[0])
    shape_data.dhdxi[3, 1] = 0.25 * (1.0 - xi[0])

    return shape_data


def get_shape_tria6(xi):
    # Check the dimension of physical space
    if len(xi) != 2:
        raise NotImplementedError('2D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(6)
    shape_data.dhdxi = empty(shape=(6, 2))
    shape_data.xi = xi

    shape_data.h[0] = 1.0 - xi[0] - xi[1]
    shape_data.h[1] = xi[0]
    shape_data.h[2] = xi[1]

    # Calculate shape functions
    shape_data.h[0] = 1.0 - xi[0] - xi[1] - 2.0 * xi[0] * (1.0 - xi[0] - xi[1]) - 2.0 * xi[1] * (1.0 - xi[0] - xi[1])
    shape_data.h[1] = xi[0] - 2.0 * xi[0] * (1.0 - xi[0] - xi[1]) - 2.0 * xi[0] * xi[1]
    shape_data.h[2] = xi[1] - 2.0 * xi[0] * xi[1] - 2.0 * xi[1] * (1.0 - xi[0] - xi[1])
    shape_data.h[3] = 4.0 * xi[0] * (1.0 - xi[0] - xi[1])
    shape_data.h[4] = 4.0 * xi[0] * xi[1]
    shape_data.h[5] = 4.0 * xi[1] * (1.0 - xi[0] - xi[1])

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -1.0 - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[0] + 2.0 * xi[1]
    shape_data.dhdxi[1, 0] = 1.0 - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[0] - 2.0 * xi[1]
    shape_data.dhdxi[2, 0] = 0.0
    shape_data.dhdxi[3, 0] = 4.0 * (1.0 - xi[0] - xi[1]) - 4.0 * xi[0]
    shape_data.dhdxi[4, 0] = 4.0 * xi[1]
    shape_data.dhdxi[5, 0] = -4.0 * xi[1]

    shape_data.dhdxi[0, 1] = -1.0 + 2.0 * xi[0] - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[1]
    shape_data.dhdxi[1, 1] = 0.0
    shape_data.dhdxi[2, 1] = 1.0 - 2.0 * xi[0] - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[1]
    shape_data.dhdxi[3, 1] = -4.0 * xi[0]
    shape_data.dhdxi[4, 1] = 4.0 * xi[0]
    shape_data.dhdxi[5, 1] = 4.0 * (1.0 - xi[0] - xi[1]) - 4.0 * xi[1]

    return shape_data


def get_shape_quad8(xi):
    # Check the dimension of physical space
    if len(xi) != 2:
        raise NotImplementedError('2D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(8)
    shape_data.dhdxi = empty(shape=(8, 2))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = -0.25 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[0] + xi[1])
    shape_data.h[1] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[0]) * (1.0 - xi[1])
    shape_data.h[2] = -0.25 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[0] + xi[1])
    shape_data.h[3] = 0.5 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[1])
    shape_data.h[4] = -0.25 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[0] - xi[1])
    shape_data.h[5] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[0]) * (1.0 + xi[1])
    shape_data.h[6] = -0.25 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[0] - xi[1])
    shape_data.h[7] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[1])

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -0.25 * (-1.0 + xi[1]) * (2.0 * xi[0] + xi[1])
    shape_data.dhdxi[1, 0] = xi[0] * (-1.0 + xi[1])
    shape_data.dhdxi[2, 0] = 0.25 * (-1.0 + xi[1]) * (-2.0 * xi[0] + xi[1])
    shape_data.dhdxi[3, 0] = -0.5 * (1.0 + xi[1]) * (-1.0 + xi[1])
    shape_data.dhdxi[4, 0] = 0.25 * (1.0 + xi[1]) * (2.0 * xi[0] + xi[1])
    shape_data.dhdxi[5, 0] = -xi[0] * (1.0 + xi[1])
    shape_data.dhdxi[6, 0] = -0.25 * (1.0 + xi[1]) * (-2.0 * xi[0] + xi[1])
    shape_data.dhdxi[7, 0] = 0.5 * (1.0 + xi[1]) * (-1.0 + xi[1])

    shape_data.dhdxi[0, 1] = -0.25 * (-1.0 + xi[0]) * (xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[1, 1] = 0.5 * (1.0 + xi[0]) * (-1.0 + xi[0])
    shape_data.dhdxi[2, 1] = 0.25 * (1.0 + xi[0]) * (-xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[3, 1] = -xi[1] * (1.0 + xi[0])
    shape_data.dhdxi[4, 1] = 0.25 * (1.0 + xi[0]) * (xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[5, 1] = -0.5 * (1.0 + xi[0]) * (-1.0 + xi[0])
    shape_data.dhdxi[6, 1] = -0.25 * (-1.0 + xi[0]) * (-xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[7, 1] = xi[1] * (-1.0 + xi[0])

    return shape_data


def get_shape_quad9(xi):
    # Check the dimension of physical space
    if len(xi) != 2:
        raise NotImplementedError('2D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(9)
    shape_data.dhdxi = empty(shape=(9, 2))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = -0.25 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[0] + xi[1])
    shape_data.h[1] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[0]) * (1.0 - xi[1])
    shape_data.h[2] = -0.25 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[0] + xi[1])
    shape_data.h[3] = 0.5 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[1])
    shape_data.h[4] = -0.25 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[0] - xi[1])
    shape_data.h[5] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[0]) * (1.0 + xi[1])
    shape_data.h[6] = -0.25 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[0] - xi[1])
    shape_data.h[7] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[1])
    shape_data.h[8] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[1])

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -0.25 * (-1.0 + xi[1]) * (2.0 * xi[0] + xi[1])
    shape_data.dhdxi[1, 0] = xi[0] * (-1.0 + xi[1])
    shape_data.dhdxi[2, 0] = 0.25 * (-1.0 + xi[1]) * (-2.0 * xi[0] + xi[1])
    shape_data.dhdxi[3, 0] = -0.5 * (1.0 + xi[1]) * (-1.0 + xi[1])
    shape_data.dhdxi[4, 0] = 0.25 * (1.0 + xi[1]) * (2.0 * xi[0] + xi[1])
    shape_data.dhdxi[5, 0] = -xi[0] * (1.0 + xi[1])
    shape_data.dhdxi[6, 0] = -0.25 * (1.0 + xi[1]) * (-2.0 * xi[0] + xi[1])
    shape_data.dhdxi[7, 0] = 0.5 * (1.0 + xi[1]) * (-1.0 + xi[1])
    shape_data.dhdxi[8, 0] = 0.5 * (1.0 + xi[1]) * (-1.0 + xi[1])

    shape_data.dhdxi[0, 1] = -0.25 * (-1.0 + xi[0]) * (xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[1, 1] = 0.5 * (1.0 + xi[0]) * (-1.0 + xi[0])
    shape_data.dhdxi[2, 1] = 0.25 * (1.0 + xi[0]) * (-xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[3, 1] = -xi[1] * (1.0 + xi[0])
    shape_data.dhdxi[4, 1] = 0.25 * (1.0 + xi[0]) * (xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[5, 1] = -0.5 * (1.0 + xi[0]) * (-1.0 + xi[0])
    shape_data.dhdxi[6, 1] = -0.25 * (-1.0 + xi[0]) * (-xi[0] + 2.0 * xi[1])
    shape_data.dhdxi[7, 1] = xi[1] * (-1.0 + xi[0])
    shape_data.dhdxi[8, 1] = xi[1] * (-1.0 + xi[0])

    return shape_data


def get_shape_tetra4(xi):
    # Check the dimension of physical space
    if len(xi) != 3:
        raise NotImplementedError('3D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(4)
    shape_data.dhdxi = empty(shape=(4, 3))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = 1.0 - xi[0] - xi[1] - xi[2]
    shape_data.h[1] = xi[0]
    shape_data.h[2] = xi[1]
    shape_data.h[3] = xi[2]

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -1.0
    shape_data.dhdxi[1, 0] = 1.0
    shape_data.dhdxi[2, 0] = 0.0
    shape_data.dhdxi[3, 0] = 0.0

    shape_data.dhdxi[0, 1] = -1.0
    shape_data.dhdxi[1, 1] = 0.0
    shape_data.dhdxi[2, 1] = 1.0
    shape_data.dhdxi[3, 1] = 0.0

    shape_data.dhdxi[0, 2] = -1.0
    shape_data.dhdxi[1, 2] = 0.0
    shape_data.dhdxi[2, 2] = 0.0
    shape_data.dhdxi[3, 2] = 1.0

    return shape_data


def get_shape_pyramid5(xi):
    # Check the dimension of physical space
    if len(xi) != 3:
        raise NotImplementedError('3D only')

    shape_data = ShapeData()

    # Set length of lists
    shape_data.h = empty(5)
    shape_data.dhdxi = empty(shape=(5, 3))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.h[1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.h[2] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.h[3] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.h[4] = 0.5 * (1.0 + xi[2])

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[1, 0] = 0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[2, 0] = 0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[3, 0] = -0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[4, 0] = 0.0

    shape_data.dhdxi[0, 1] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[1, 1] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[2, 1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[3, 1] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[4, 1] = 0.0

    shape_data.dhdxi[0, 2] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[1])
    shape_data.dhdxi[1, 2] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[1])
    shape_data.dhdxi[2, 2] = -0.125 * (1.0 + xi[0]) * (1.0 + xi[1])
    shape_data.dhdxi[3, 2] = -0.125 * (1.0 - xi[0]) * (1.0 + xi[1])
    shape_data.dhdxi[4, 2] = 0.5

    return shape_data


def get_shape_prism6(xi):
    # Check the dimension of physical space
    if len(xi) != 3:
        raise NotImplementedError('3D only')

    # Initialise tuples

    shape_data = ShapeData()

    shape_data.h = empty(6)
    shape_data.dhdxi = empty(shape=(6, 3))
    shape_data.xi = xi

    shape_data_line2 = get_shape_line2(xi[2])
    shape_data_tria3 = get_shape_tria3(xi[:2])

    for i in range(3):
        for j in range(2):
            shape_data.h[i * 2 + j] = shape_data_line2.h[j] * shape_data_tria3.h[i]

            shape_data.dhdxi[i * 2 + j, 0] = shape_data_line2.h[j] * shape_data_tria3.dhdxi[i, 0]
            shape_data.dhdxi[i * 2 + j, 1] = shape_data_line2.h[j] * shape_data_tria3.dhdxi[i, 1]
            shape_data.dhdxi[i * 2 + j, 2] = shape_data_line2.dhdxi[j, 0] * shape_data_tria3.h[i]

    return shape_data


def get_shape_prism18(xi):
    # Check the dimension of physical space
    if len(xi) != 3:
        raise NotImplementedError('3D only')

    shape_data = ShapeData()
    shape_data.h = empty(18)
    shape_data.dhdxi = empty(shape=(18, 3))
    shape_data.xi = xi

    shape_data_line3 = get_shape_line3(xi[3])
    shape_data_tria6 = get_shape_tria6(xi[:2])

    for i in range(6):
        for j in range(3):
            shape_data.h[i * 3 + j] = shape_data_line3.h[j] * shape_data_tria6.h[i]

            shape_data.dhdxi[i * 3 + j, 0] = shape_data_line3.h[j] * shape_data_tria6.dhdxi[i, 0]
            shape_data.dhdxi[i * 3 + j, 1] = shape_data_line3.h[j] * shape_data_tria6.dhdxi[i, 1]
            shape_data.dhdxi[i * 3 + j, 2] = shape_data_line3.dhdxi[j, 0] * shape_data_tria6.h[i]

    return shape_data


def get_shape_hex8(xi):
    if len(xi) != 3:
        raise NotImplementedError('The isoparamatric coordinate should be 3D.')

    shape_data = ShapeData()

    shape_data.h = empty(8)
    shape_data.dhdxi = empty(shape=(8, 3))
    shape_data.xi = xi

    # Calculate shape functions
    shape_data.h[0] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.h[1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.h[2] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.h[3] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.h[4] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
    shape_data.h[5] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
    shape_data.h[6] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])
    shape_data.h[7] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])

    # Calculate derivatives of shape functions
    shape_data.dhdxi[0, 0] = -0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[1, 0] = 0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[2, 0] = 0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[3, 0] = -0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    shape_data.dhdxi[4, 0] = -0.125 * (1.0 - xi[1]) * (1.0 + xi[2])
    shape_data.dhdxi[5, 0] = 0.125 * (1.0 - xi[1]) * (1.0 + xi[2])
    shape_data.dhdxi[6, 0] = 0.125 * (1.0 + xi[1]) * (1.0 + xi[2])
    shape_data.dhdxi[7, 0] = -0.125 * (1.0 + xi[1]) * (1.0 + xi[2])

    shape_data.dhdxi[0, 1] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[1, 1] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[2, 1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[3, 1] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    shape_data.dhdxi[4, 1] = -0.125 * (1.0 - xi[0]) * (1.0 + xi[2])
    shape_data.dhdxi[5, 1] = -0.125 * (1.0 + xi[0]) * (1.0 + xi[2])
    shape_data.dhdxi[6, 1] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[2])
    shape_data.dhdxi[7, 1] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[2])

    shape_data.dhdxi[0, 2] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[1])
    shape_data.dhdxi[1, 2] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[1])
    shape_data.dhdxi[2, 2] = -0.125 * (1.0 + xi[0]) * (1.0 + xi[1])
    shape_data.dhdxi[3, 2] = -0.125 * (1.0 - xi[0]) * (1.0 + xi[1])
    shape_data.dhdxi[4, 2] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1])
    shape_data.dhdxi[5, 2] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1])
    shape_data.dhdxi[6, 2] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1])
    shape_data.dhdxi[7, 2] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1])

    return shape_data


def get_element_type(element_coords):
    num_element_nodes = element_coords.shape[0]
    rank = element_coords.shape[1]

    if rank == 1:
        if num_element_nodes == 2:
            return "line2"
        elif num_element_nodes == 3:
            return "line3"
        else:
            raise NotImplementedError('No 1D element with ' + str(num_element_nodes) + ' nodes available')
    elif rank == 2:
        if num_element_nodes == 3:
            return "tria3"
        elif num_element_nodes == 4:
            return "quad4"
        elif num_element_nodes == 6:
            return "tria6"
        elif num_element_nodes == 8:
            return "quad8"
        elif num_element_nodes == 9:
            return "quad9"
        else:
            raise NotImplementedError('No 2D element with ' + str(num_element_nodes) + ' nodes available')
    elif rank == 3:
        if num_element_nodes == 4:
            return "tetra4"
        elif num_element_nodes == 5:
            return "pyramid5"
        elif num_element_nodes == 6:
            return "prism6"
        elif num_element_nodes == 8:
            return "hex8"
        elif num_element_nodes == 18:
            return "prism18"
        else:
            raise NotImplementedError('No 3D element with ' + str(num_element_nodes) + ' nodes available')
    else:
        raise NotImplementedError('Rank must be 1, 2 or 3')


def tria_scheme(order):
    if order == 1:
        xi = [[1.0 / 3.0, 1.0 / 3.0]]
        weight = [0.5]
    elif order == 3:
        r1 = 1.0 / 6.0
        r2 = 2.0 / 3.0

        xi = [[r1, r1], [r2, r1], [r1, r2]]

        w1 = 1.0 / 6.0

        weight = [w1, w1, w1]
    elif order == 7:
        r1 = 0.5 * 0.1012865073235
        r2 = 0.5 * 0.7974269853531
        r4 = 0.5 * 0.4701420641051
        r6 = 0.0597158717898
        r7 = 1.0 / 3.0

        xi = [[r1, r1], [r2, r1], [r1, r2], [r4, r6], [r4, r4], [r6, r4], [r7, r7]]

        w1 = 0.1259391805448
        w4 = 0.1323941527885
        w7 = 0.225

        weight = [w1, w1, w1, w4, w4, w4, w7]
    else:
        raise NotImplementedError('Order must be 1, 3 or 7')

    return xi, weight


def tetra_scheme(order):
    if order == 1:
        third = 1. / 3.

        xi = [[third, third, third]]
        weight = [0.5 * third]
    else:
        raise NotImplementedError('Only order 1 integration implemented')

    return xi, weight


def pyramid_scheme(order):
    if order == 1:
        xi = [[0., 0., -0.5]]
        weight = [128.0 / 27.0]  # [8.0/3.0] #[ 18.967 ]
    else:
        raise NotImplementedError('Only order 1 integration implemented')

    return xi, weight


def get_integration_points(element_type, order, method):
    xi = []
    weight = []

    if element_type[:-1] == "line":
        if element_type == "line2":
            standard_order = 2
        elif element_type == "line3":
            standard_order = 3
        else:
            raise NotImplementedError('Unsupported ' + element_type)
        xi, weight = gauss_scheme(standard_order + order)
        xi = [float(a.real) for a in xi]

    elif element_type[:-1] == "tria":
        order_array = [1, 3, 7]
        if element_type == "tria3":
            standard_order = 0
        elif element_type == "tria6":
            standard_order = 1
        else:
            raise NotImplementedError('Unsupported ' + element_type)
        xi, weight = tria_scheme(order_array[standard_order + order])

    elif element_type[:-1] == "tetra":
        standard_order = 1
        xi, weight = tetra_scheme(standard_order + order)

    elif element_type == "pyramid5":
        standard_order = 1
        xi, weight = pyramid_scheme(standard_order + order)

    elif element_type[:-1] == "quad":
        if element_type == "quad4":
            standard_order = 2
        elif element_type == "quad8" or element_type == "quad9":
            standard_order = 3
        else:
            raise NotImplementedError('Unsupported ' + element_type)
        standard_order += order

        ip, w = gauss_scheme(standard_order)

        for i in range(standard_order):
            for j in range(standard_order):
                xi.append([float(ip[i].real), float(ip[j].real)])
                weight.append(w[i] * w[j])

    elif element_type[:-1] == "hex":
        if element_type == "hex8":
            standard_order = 2
        else:
            raise NotImplementedError('Unsupported ' + element_type)

        standard_order += order

        ip, w = gauss_scheme(standard_order)

        for i in range(standard_order):
            for j in range(standard_order):
                for k in range(standard_order):
                    xi.append([float(ip[i].real), float(ip[j].real), float(ip[k].real)])
                    weight.append(w[i] * w[j] * w[k])

    elif element_type[:-1] == "prism":
        order_array = [1, 3, 7]
        if element_type == "prism6":
            standard_order = 2
        elif element_type == "prism18":
            standard_order = 3
        else:
            raise NotImplementedError('Unsupported ' + element_type)

        standard_order += order

        ip0, w0 = tria_scheme(order_array[standard_order])
        ip1, w1 = gauss_scheme(standard_order)

        for i in range(order_array[standard_order]):
            for j in range(standard_order):
                xi.append([float(ip0[i][0].real), float(ip0[i][1].real), float(ip1[j].real)])
                weight.append(w0[i] * w1[j])

    return xi, weight


def calc_weight_and_derivatives(element_coords, shape_data, weight):
    jac = dot(element_coords.transpose(), shape_data.dhdxi)

    if jac.shape[0] == jac.shape[1]:
        shape_data.dhdx = dot(shape_data.dhdxi, inv(jac))
        shape_data.weight = abs(det(jac)) * weight

    elif jac.shape[0] == 2 and jac.shape[1] == 1:
        shape_data.weight = sqrt(sum(sum(jac * jac))) * weight

    elif jac.shape[0] == 3 and jac.shape[1] == 2:
        jac3 = zeros(shape=(3, 3))

        jac3[:, :2] = jac

        dA = zeros(3)

        dA[0] = norm(cross(jac3[:, 1], jac3[:, 2]))
        dA[1] = norm(cross(jac3[:, 2], jac3[:, 0]))
        dA[2] = norm(cross(jac3[:, 0], jac3[:, 1]))

        shape_data.weight = norm(dA) * weight


def get_element_shape_data(element_coords, order=0, method='Gauss', element_type='Default'):
    element_data = ElementShapeData()

    if element_type == 'Default':
        element_type = get_element_type(element_coords)

    (ip_coords, ip_wights) = get_integration_points(element_type, order, method)

    for xi, weight in zip(ip_coords, ip_wights):
        try:
            shape_data = eval('get_shape_' + element_type + '(xi)')
        except NotImplementedError:
            raise NotImplementedError('Unknown type :' + element_type)

        calc_weight_and_derivatives(element_coords, shape_data, weight)

        shape_data.x = dot(shape_data.h, element_coords)

        element_data.shape_data.append(shape_data)

    return element_data


def get_shape_data(order=0, method='Gauss', element_type='Default'):
    shape_data = ElementShapeData()

    (ip_coords, ip_wights) = get_integration_points(element_type, order, method)

    for xi, weight in zip(ip_coords, ip_wights):
        try:
            shape_data = eval('getShape' + element_type + '(xi)')
        except NotImplementedError:
            raise NotImplementedError('Unknown type :' + element_type)

        shape_data.dhdx = shape_data.dhdxi
        shape_data.weight = weight

        shape_data.shape_data.append(shape_data)

    return shape_data
