from numpy import dot, empty, real, sqrt
from scipy.linalg import det, inv

from .shape_functions import ShapeData, ElementShapeData, get_integration_points, get_element_type




def getBezierLine4(xi, C):
    # Check the dimensions of the parametric space

    if type(xi) == 'numpy.float64':
        raise NotImplementedError('1D only')
    if C.shape[1] != 4:
        raise NotImplementedError('C needs to have 4 columns.')

    shape_data = ShapeData()

    # Set length of lists

    #  shape_data.h     = empty( 4 )
    #  shape_data.dhdxi = empty( shape=(1,4) )
    shape_data.xi = xi

    B = empty(4)
    dBdxi = empty(shape=(4, 1))

    # Calculate shape functions

    B[0] = -0.125 * (xi - 1.) ** 3
    B[1] = 0.375 * (xi - 1.) ** 2 * (xi + 1.)
    B[2] = -0.375 * (xi - 1.) * (xi + 1.) ** 2
    B[3] = 0.125 * (xi + 1.) ** 3

    # Calculate derivatives of shape functions

    dBdxi[0, 0] = -0.375 * (xi - 1.) ** 2
    dBdxi[1, 0] = 0.75 * (xi - 1.0) * (xi + 1.0) + 0.375 * (xi - 1.) ** 2
    dBdxi[2, 0] = -0.375 * (1 + xi) ** 2 - 0.75 * (1 + xi) * (xi - 1)
    dBdxi[3, 0] = 0.375 * (xi + 1.) ** 2

    shape_data.h = dot(C, B)
    shape_data.dhdxi = dot(C, dBdxi)

    return shape_data


# -------------

def calcWeight(jac):
    n = jac.shape

    if n[0] == n[1]:
        return det(jac)
    elif n[0] == 1 and n[1] == 2:
        return sqrt(sum(sum(jac * jac)))




def getElemBezierData(element_coords, C, order=4, method="Gauss", element_type='default'):
    element_data = ElementShapeData()

    if element_type == 'default':
        element_type = get_element_type(element_coords)

    (gp_coords, gp_wights) = get_integration_points("Line3", order, method)

    for xi, intWeight in zip(real(gp_coords), gp_wights):
        try:
            shape_data = eval('getBezier' + element_type + '(xi,C)')
        except:
            raise NotImplementedError('Unknown type :' + element_type)

        jac = dot(shape_data.dhdxi.transpose(), element_coords)

        if jac.shape[0] is jac.shape[1]:
            shape_data.dhdx = (dot(inv(jac), shape_data.dhdxi.transpose())).transpose()

        shape_data.weight = calcWeight(jac) * intWeight

        element_data.shape_data.append(shape_data)

    return element_data
