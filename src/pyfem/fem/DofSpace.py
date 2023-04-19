from copy import deepcopy

import scipy.linalg
from numpy import array, dot, zeros, where
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import spsolve

from pyfem.fem.Constrainer import Constrainer
from pyfem.utils.parser import read_node_table
from pyfem.utils.itemList import itemList
from pyfem.utils.logger import getLogger

logger = getLogger()


class DofSpace:
    '''
    Class dofspace
    '''

    def __init__(self, elements):

        '''
        Constructor
        '''

        self.dofTypes = elements.getDofTypes()
        self.dofs = array(list(range(len(elements.nodes) * len(self.dofTypes)))).reshape(
            (len(elements.nodes), len(self.dofTypes)))
        self.nodes = elements.nodes

        # Create the ID map
        self.IDmap = itemList()
        for ind, ID in enumerate(elements.nodes):
            self.IDmap.add(ID, ind)

        self.allConstrainedDofs = []



    def __str__(self):

        '''
        Prints the total overview of degrees of freedom
        '''

        return str(self.dofs)



    def __len__(self):

        '''
        Function that returns the length of the dofpsace, i.e. the number of
        degrees of freedeom
        '''

        return len(self.dofs.flatten())



    def setConstrainFactor(self, fac, loadCase="All_"):

        if loadCase == "All_":
            for name in self.cons.constrainedFac.keys():
                self.cons.constrainedFac[name] = fac
        else:
            self.cons.constrainedFac[loadCase] = fac



    def readFromFile(self, fname):

        logger.info("Reading constraints ..........")

        NodeTable = read_node_table(fname, "NodeConstraints", self.nodes)

        self.cons = self.createConstrainer(NodeTable)



    def createConstrainer(self, nodeTables=None):

        cons = Constrainer(len(self))

        if nodeTables == None:
            label = "main"
            cons.constrainedDofs[label] = []
            cons.constrainedVals[label] = []
            cons.constrainedFac[label] = 1.0

            self.cons = cons
            return cons

        for NodeTable in nodeTables:

            label = NodeTable.sub_label

            cons.constrainedDofs[label] = []
            cons.constrainedVals[label] = []
            cons.constrainedFac[label] = 1.0

            for item in NodeTable.data:

                node_id = item[1]
                dof_type = item[0]
                val = item[2]

                if not node_id in self.nodes:
                    raise RuntimeError('Node ID ' + str(node_id) + ' does not exist')

                ind = self.IDmap.get(node_id)

                if dof_type not in self.dofTypes:
                    raise RuntimeError('DOF type "' + dof_type + '" does not exist')

                if len(item) == 3:
                    dofID = self.dofs[ind, self.dofTypes.index(dof_type)]

                    cons.addConstraint(dofID, val, label)
                else:
                    slaveNodeID = item[4]
                    slaveDofType = item[3]
                    factor = item[5]

                    if not slaveNodeID[0] in self.nodes:
                        raise RuntimeError('Node ID ' + str(slaveNodeID) + ' does not exist')

                    slaveInd = self.IDmap.get(slaveNodeID)

                    if slaveDofType not in self.dofTypes:
                        raise RuntimeError('DOF type "' + slaveDofType + '" does not exist')

                    slaveDof = self.dofs[slaveInd, self.dofTypes.index(slaveDofType)]

                    dofID = self.dofs[ind, self.dofTypes.index(dof_type)]

                    cons.addConstraint(dofID, [val, slaveDof, factor], label)

                    # Check for all tyings whether master of slave is not slave itself
        cons.checkConstraints(self, nodeTables)

        cons.flush()

        return cons



    def getForType(self, node_ids, dof_type):

        '''
        Returns all dofIDs for given dof_type for a list of nodes
        '''

        return self.dofs[self.IDmap.get(node_ids), self.dofTypes.index(dof_type)]



    def getForTypes(self, node_ids, dofTypes):

        '''
        Returns all dofIDs for given list of dof_type for a list of nodes
        '''

        dofs = []

        for node in node_ids:
            for dof_type in dofTypes:
                dofs.append(self.dofs[self.IDmap.get(node), self.dofTypes.index(dof_type)])

        return dofs



    def getDofName(self, dofID):

        '''
        Returns the dofID as a string. For example 'u[14]'
        '''

        return self.getTypeName(dofID) + '[' + str(self.getNodeID(dofID)) + ']'



    def getNodeID(self, dofID):

        '''
        Returns the node ID of dofID
        '''

        return self.nodes.findID(int(where(self.dofs == dofID)[0]))



    def get_type(self, dofID):

        '''
        Returns the type of dofID
        '''

        return int(where(self.dofs == dofID)[1])



    def getTypeName(self, dofID):

        '''
        Returns the name of the dof_type
        '''

        return self.dofTypes[self.get_type(dofID)]



    def get(self, node_ids):

        '''Returns all dofIDs for a list of nodes'''

        return self.dofs[self.IDmap.get(node_ids)].flatten()



    def copyConstrainer(self, dofTypes: list = None):

        '''
        
        '''

        newCons = deepcopy(self.cons)

        if type(dofTypes) is str:
            dofTypes = [dofTypes]

        for dof_type in dofTypes:
            for iDof in self.dofs[:, self.dofTypes.index(dof_type)]:
                for label in newCons.constrainedFac.keys():
                    newCons.addConstraint(iDof, 0.0, label)

        newCons.flush()

        return newCons

    # -------------------------------------------------------------------------------
    #  
    # -------------------------------------------------------------------------------

    def solve(self, A, b, constrainer=None):

        '''Solves the system Ax = b using the internal constraints matrix.
           Returns the total solution vector x.'''

        if constrainer is None:
            constrainer = self.cons

        if len(A.shape) == 2:

            a = zeros(len(self))

            constrainer.addConstrainedValues(a)

            A_constrained = constrainer.C.transpose() * (A * constrainer.C)

            b_constrained = constrainer.C.transpose() * (b - A * a)

            x_constrained = spsolve(A_constrained, b_constrained)

            x = constrainer.C * x_constrained

            constrainer.addConstrainedValues(x)

        elif len(A.shape) == 1:
            x = b / A

            constrainer.setConstrainedValues(x)

        return x



    def eigensolve(self, A, B, count=5):

        '''Calculates the first count eigenvalues and eigenvectors of a
           system with ( A lambda B ) x '''

        A_constrained = dot(dot(self.cons.C.transpose(), A), self.cons.C)
        B_constrained = dot(dot(self.cons.C.transpose(), B), self.cons.C)

        eigvals, eigvecs = eigsh(A_constrained, count, B_constrained, sigma=0., which='LM')

        x = zeros(shape=(len(self), count))

        for i, psi in enumerate(eigvecs.transpose()):
            x[:, i] = self.cons.C * psi

        return eigvals, x



    def norm(self, r, constrainer=None):

        '''
        Calculates the norm of vector r excluding the constrained dofs
        '''

        if constrainer is None:
            constrainer = self.cons

        return scipy.linalg.norm(constrainer.C.transpose() * r)



    def maskPrescribed(self, a, val=0.0, constrainer=None):

        '''
        Replaced the prescribed dofs by val
        '''

        if constrainer is None:
            constrainer = self.cons

        a[constrainer.constrainedDofs["None"]] = val

        return a
