from numpy import zeros

from pyfem.materials.MaterialManager import MaterialManager




class Element(list):
    dofTypes = []



    def __init__(self, elnodes, props):
        list.__init__(self, elnodes)

        self.family = "CONTINUUM"
        self.history = {}
        self.current = {}
        self.solverStat = props.solverStat

        for name, val in props:
            if name == "material":
                self.matProps = val

                self.matProps.rank = props.rank
                self.matProps.solverStat = self.solverStat
                self.mat = MaterialManager(self.matProps)

            setattr(self, name, val)



    def dofCount(self):

        return len(self) * len(self.dofTypes)



    def getNodes(self):
        return self



    def get_type(self):
        return self.element_type



    def appendNodalOutput(self, labels, data, weight=1.0):

        for i, name in enumerate(labels):
            if not hasattr(self.globdat, name):
                self.globdat.outputNames.append(name)

                setattr(self.globdat, name, zeros(len(self.globdat.nodes)))
                setattr(self.globdat, name + 'Weights', zeros(len(self.globdat.nodes)))

            outMat = getattr(self.globdat, name)
            outWeights = getattr(self.globdat, name + 'Weights')

            if data.ndim == 1:
                for idx in self.globdat.nodes.getIndices(self):
                    outMat[idx] += data[i]
                    outWeights[idx] += weight
            else:
                for j, idx in enumerate(self.globdat.nodes.getIndices(self)):
                    outMat[idx] += data[j, i]
                    outWeights[idx] += weight



    def setHistoryParameter(self, name, val):
        self.current[name] = val



    def getHistoryParameter(self, name):
        return self.history[name]



    def commitHistory(self):
        self.history = self.current.copy()
        self.current = {}

        if hasattr(self, "mat"):
            self.mat.commitHistory()

    def commit(self, elemdat):
        pass

    #
    #
    #

    def loadFactor(self):

        return self.solverStat.lam
