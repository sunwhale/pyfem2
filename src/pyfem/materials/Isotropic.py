from numpy import zeros, dot

from pyfem.materials.BaseMaterial import BaseMaterial


class Isotropic(BaseMaterial):

    def __init__(self, props):

        self.incremental = False

        # Call the BaseMaterial constructor
        BaseMaterial.__init__(self, props)

        # Create the hookean matrix
        self.H = zeros((6, 6))

        fac = 1.0 / (2.0 * self.nu * self.nu + self.nu - 1.0);

        self.H[0, 0] = fac * self.E * (self.nu - 1.0);
        self.H[0, 1] = -1.0 * fac * self.E * self.nu;
        self.H[0, 2] = self.H[0, 1];
        self.H[1, 0] = self.H[0, 1];
        self.H[1, 1] = self.H[0, 0];
        self.H[1, 2] = self.H[0, 1];
        self.H[2, 0] = self.H[0, 1];
        self.H[2, 1] = self.H[0, 1];
        self.H[2, 2] = self.H[0, 0];
        self.H[3, 3] = self.E / (2.0 + 2.0 * self.nu);
        self.H[4, 4] = self.H[3, 3];
        self.H[5, 5] = self.H[3, 3];

        # Set the labels for the output data in this material model
        self.outLabels = ["S11", "S22", "S33", "S23", "S13", "S12"]

        if self.incremental:
            self.setHistoryParameter('sigma', zeros(6))
            self.commit_history()

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def getStress(self, deformation):

        if self.incremental:
            sigma = self.getHistoryParameter('sigma')
            sigma += dot(self.H, deformation.dstrain)
            self.setHistoryParameter('sigma', sigma)
        else:
            sigma = dot(self.H, deformation.strain)

        self.outData = sigma

        return sigma, self.H

    def getTangent(self):

        return self.H
