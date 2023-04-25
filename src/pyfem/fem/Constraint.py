from numpy import array
from scipy.sparse import coo_matrix

from pyfem.utils.logger import get_logger

logger = get_logger()


class Constraint:

    def __init__(self, nDofs, name="Main"):

        self.nDofs = nDofs
        self.constrainData = {}
        self.name = name

        self.constrained_dofs = {}
        self.constrained_values = {}
        self.constrained_factors = {}

    def add_constraint(self, dof_id, val, label):

        if dof_id in self.constrainData:
            self.constrainData[dof_id].append(val)
        else:
            self.constrainData[dof_id] = [val]

        if (type(val) is list) and (len(val) == 3):
            addVal = val[0]
        else:
            addVal = val

        if dof_id in self.constrained_dofs[label]:
            self.setFactorForDof(addVal, dof_id, label)
            return

        self.constrained_dofs[label].append(dof_id)

        self.constrained_values[label].append(addVal)

        logger.debug("TEST")

    def check_constraints(self, dofspace, node_tables):

        '''Checks tying relations between dofs'''
        for item in self.constrainData:

            dofInd = item

            for ilabel in self.constrained_dofs:
                if dofInd in self.constrained_dofs[ilabel]:
                    label = ilabel

            removeItem = []

            for tiedItem in self.constrainData[item]:

                if (type(tiedItem) is list) and len(tiedItem) == 3:
                    Dat = tiedItem

                    valSlave = Dat[0]
                    masterDofID = Dat[1][0]
                    facSlave = Dat[2]

                    if masterDofID in self.constrainData:
                        tempVal = []
                        tempFac = []

                        # Recursive loop until masterDofID not a list, but prescribed value         
                        while masterDofID in self.constrainData:
                            master = self.constrainData[masterDofID][0]
                            if type(master) is list and len(master) == 3:
                                masterDofID = master[1][0]
                                tempVal.append(master[0])
                                tempFac.append(master[2])
                            else:
                                masterFin = master
                                masterDofID = -1

                        for iVal, iFac in reversed(list(zip(tempVal, tempFac))):
                            masterFin += iVal + master * iFac

                        self.add_constraint(dofInd, valSlave + masterFin * facSlave, label)

                        removeItem.append(tiedItem)

            for iRemove in removeItem:
                self.constrainData[item].remove(iRemove)

    def flush(self):

        '''Returns the constraints matrix using the class member arrays'''

        row = []
        col = []
        val = []
        master = {}

        iCon = 0

        for dof_id in range(self.nDofs):
            if dof_id in self.constrainData:
                for item in self.constrainData[dof_id]:
                    if type(item) is not list:
                        continue
                    else:
                        if item[1][0] in self.constrainData:
                            # Something not checked correct in checkConstraint2
                            raise RuntimeError('ERROR - Master of slave is a slave itself')
                        else:
                            master[dof_id] = item

            else:
                row.append(dof_id)
                col.append(iCon)
                val.append(1.)

                iCon += 1

        # Assign correct slaves to masters
        for iSlave in master:
            fac = master[iSlave][2]
            masterDofID = master[iSlave][1][0]

            # Find column for free DOF of the Master
            rowID = row.index(masterDofID)
            col.append(rowID)
            row.append(iSlave)
            val.append(fac)

        self.C = coo_matrix((val, (row, col)), shape=(self.nDofs, iCon))

    def add_constrained_values(self, a):

        for name in self.constrained_dofs.keys():
            a[self.constrained_dofs[name]] += self.constrained_factors[name] * array(self.constrained_values[name])

    def set_constrained_values(self, a):

        for name in self.constrained_dofs.keys():
            a[self.constrained_dofs[name]] = self.constrained_factors[name] * array(self.constrained_values[name])

    def set_constrain_factor(self, fac, load_case="All_"):

        if load_case == "All_":
            for name in self.constrained_factors.keys():
                self.constrained_factors[name] = fac
        else:
            self.constrained_factors[load_case] = fac

    def setPrescribedDofs(self, a, val=0.0):

        for name in self.constrained_dofs.keys():
            a[self.constrained_dofs[name]] = array(val)

    def setFactorForDof(self, fac, dof_id, label):

        idx = self.constrained_dofs[label].index(dof_id)
        self.constrained_values[label][idx] += fac

    def slaveCount(self):

        counter = 0
        for name in self.constrained_factors.keys():
            counter += len(self.constrained_dofs[name])

        return counter
