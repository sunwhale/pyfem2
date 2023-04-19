from numpy import zeros, dot

from pyfem.utils.kinematics import Kinematics
from pyfem.utils.shape_functions import get_element_shape_data
from .Element import Element


class SmallStrainContinuum(Element):

    def __init__(self, elnodes, props):
        Element.__init__(self, elnodes, props)

        self.rank = props.rank

        if self.rank == 2:
            self.dofTypes = ['u', 'v']
            self.nstr = 3
            self.outputLabels = ["s11", "s22", "s12"]
        elif self.rank == 3:
            self.dofTypes = ['u', 'v', 'w']
            self.nstr = 6
            self.outputLabels = ["s11", "s22", "s33", "s23", "s13", "s12"]

        self.kin = Kinematics(self.rank, self.nstr)

    def __type__(self):
        return name

    # ------------------------------------------------------------------------

    def getTangentStiffness(self, elemdat):

        shape_data = get_element_shape_data(elemdat.coords)

        elemdat.outlabel.append(self.outputLabels)
        elemdat.outdata = zeros(shape=(len(elemdat.nodes), self.nstr))

        for iData in shape_data:
            b = self.getBmatrix(iData.dhdx)

            self.kin.strain = dot(b, elemdat.state)
            self.kin.dstrain = dot(b, elemdat.Dstate)

            sigma, tang = self.mat.getStress(self.kin)

            elemdat.stiff += dot(b.transpose(), dot(tang, b)) * iData.weight
            elemdat.fint += dot(b.transpose(), sigma) * iData.weight

            self.appendNodalOutput(self.mat.outLabels(), self.mat.outData())

    # -------------------------------------------------------------------------

    def getInternalForce(self, elemdat):

        shape_data = get_element_shape_data(elemdat.coords)

        elemdat.outlabel.append(self.outputLabels)
        elemdat.outdata = zeros(shape=(len(elemdat.nodes), self.nstr))

        for iData in shape_data:
            b = self.getBmatrix(iData.dhdx)

            self.kin.strain = dot(b, elemdat.state)
            self.kin.dstrain = dot(b, elemdat.Dstate)

            sigma, tang = self.mat.getStress(self.kin)

            elemdat.fint += dot(b.transpose(), sigma) * iData.weight

            self.appendNodalOutput(self.mat.outLabels(), self.mat.outData())

    # -------------------------------------------------------------------------------

    def getDissipation(self, elemdat):

        shape_data = get_element_shape_data(elemdat.coords)

        for iData in shape_data:
            b = self.getBmatrix(iData.dhdx)

            self.kin.strain = dot(b, elemdat.state)
            self.kin.dstrain = dot(b, elemdat.Dstate)

            self.mat.getStress(self.kin)

            self.kin.dgdstrain = zeros(3)
            self.kin.g = 0.0

            elemdat.fint += dot(b.transpose(), self.kin.dgdstrain) * iData.weight
            elemdat.diss += self.kin.g * iData.weight

    

    def getMassMatrix(self, elemdat):

        shape_data = get_element_shape_data(elemdat.coords)

        rho = elemdat.matprops.rho

        for iData in shape_data:
            N = self.getNmatrix(iData.h)
            elemdat.mass += dot(N.transpose(), N) * rho * iData.weight

        elemdat.lumped = sum(elemdat.mass)

    # --------------------------------------------------------------------------

    def getBmatrix(self, dphi):

        b = zeros(shape=(self.nstr, self.dofCount()))

        if self.rank == 2:
            for i, dp in enumerate(dphi):
                b[0, i * 2 + 0] = dp[0]
                b[1, i * 2 + 1] = dp[1]
                b[2, i * 2 + 0] = dp[1]
                b[2, i * 2 + 1] = dp[0]
        elif self.rank == 3:
            for i, dp in enumerate(dphi):
                b[0, i * 3 + 0] = dp[0]
                b[1, i * 3 + 1] = dp[1]
                b[2, i * 3 + 2] = dp[2]

                b[3, i * 3 + 1] = dp[2]
                b[3, i * 3 + 2] = dp[1]

                b[4, i * 3 + 0] = dp[2]
                b[4, i * 3 + 2] = dp[0]

                b[5, i * 3 + 0] = dp[1]
                b[5, i * 3 + 1] = dp[0]

        return b

    # ------------------------------------------------------------------------------

    def getNmatrix(self, h):

        N = zeros(shape=(self.rank, self.rank * len(h)))

        for i, a in enumerate(h):
            for j in list(range(self.rank)):
                N[j, self.rank * i + j] = a

        return N
