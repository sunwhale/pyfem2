from numpy import zeros

from pyfem.utils.logger import getLogger

logger = getLogger()


def clean_variable(a):
    if a == 'true':
        return True
    elif a == 'false':
        return False
    else:
        try:
            return eval(a)
        except:
            return a


class SolverStatus:

    def __init__(self):
        self.cycle = 0
        self.iiter = 0
        self.time = 0.0
        self.time0 = 0.0
        self.dtime = 0.0
        self.lam = 1.0

    def increaseStep(self):
        self.cycle += 1
        self.time += self.dtime
        self.iiter = 0


class Properties:

    def __init__(self, dictionary=None):

        if dictionary is None:
            dictionary = {}
        for key in dictionary.keys():
            setattr(self, key, dictionary[key])

    def __str__(self):

        msg = ''
        for att in dir(self):

            # Ignore private members and standard routines
            if att.startswith('__'):
                continue

            msg += 'Attribute: ' + att + '\n'
            msg += str(getattr(self, att)) + '\n'

        return msg

    def __iter__(self):

        props_list = []
        for att in dir(self):

            # Ignore private members and standard routines
            if att.startswith('__'):
                continue

            props_list.append((att, getattr(self, att)))

        return iter(props_list)

    def store(self, key, val):

        if not '.' in key:
            setattr(self, key, val)
        else:
            kets = key.split(".")

            props = self
            for y in kets[:-1]:
                props = getattr(props, y)

            setattr(props, kets[-1], clean_variable(val))


class GlobalData(Properties):

    def __init__(self, nodes, elements, dofs):

        Properties.__init__(self, {'nodes': nodes, 'elements': elements, 'dofs': dofs})

        self.state = zeros(len(self.dofs))
        self.Dstate = zeros(len(self.dofs))
        self.fint = zeros(len(self.dofs))
        self.fhat = zeros(len(self.dofs))

        self.velo = zeros(len(self.dofs))
        self.acce = zeros(len(self.dofs))

        self.SolverStatus = elements.solverStat

        self.outputNames = []

    def readFromFile(self, fname):

        logger.info("Reading external forces ......")

        fin = open(fname)

        while True:
            line = fin.readline()

            if line.startswith('<ExternalForces>') == True:
                while True:
                    line = fin.readline()

                    if line.startswith('</ExternalForces>') == True:
                        return

                    a = line.strip().split(';')

                    if len(a) == 2:
                        b = a[0].split('=')

                        if len(b) == 2:
                            c = b[0].split('[')

                            dof_type = c[0]
                            node_id = eval(c[1].split(']')[0])

                            self.fhat[self.dofs.getForType(node_id, dof_type)] = eval(b[1])

    def printNodes(self, file_name=None, inodes=None):

        if file_name is None:
            f = None
        else:
            f = open(file_name, "w")

        if inodes is None:
            inodes = list(self.nodes.keys())

        print('   Node | ', file=f, end=' ')

        for dof_type in self.dofs.dofTypes:
            print("  %-10s" % dof_type, file=f, end=' ')

        if hasattr(self, 'fint'):
            for dof_type in self.dofs.dofTypes:
                print(" fint-%-6s" % dof_type, file=f, end=' ')

        for name in self.outputNames:
            print(" %-11s" % name, file=f, end=' ')

        print(" ", file=f)
        print(('-' * 100), file=f)

        for node_id in inodes:
            print('  %4i  | ' % node_id, file=f, end=' ')
            for dof_type in self.dofs.dofTypes:
                print(' %10.3e ' % self.state[self.dofs.getForType(node_id, dof_type)], file=f, end=' ')
            for dof_type in self.dofs.dofTypes:
                print(' %10.3e ' % self.fint[self.dofs.getForType(node_id, dof_type)], file=f, end=' ')

            for name in self.outputNames:
                print(' %10.3e ' % self.getData(name, node_id), file=f, end=' ')

            print(" ", file=f)
        print(" ", file=f)

        if file_name is not None:
            f.close()

    def getData(self, outputName, inodes):

        data = getattr(self, outputName)
        weights = getattr(self, outputName + 'Weights')

        if type(inodes) is int:
            i = list(self.nodes.keys()).index(inodes)
            return data[i] / weights[i]
        else:
            outdata = []

            for row, w in zip(data[inodes], weights[inodes]):
                if w != 0:
                    outdata.append(row / w)
                else:
                    outdata.append(row)

            return outdata

    def resetNodalOutput(self):

        for outputName in self.outputNames:
            delattr(self, outputName)
            delattr(self, outputName + 'Weights')

        self.outputNames = []


class elementData():

    def __init__(self, elstate, elDstate):
        nDof = len(elstate)

        self.state = elstate
        self.Dstate = elDstate
        self.stiff = zeros(shape=(nDof, nDof))
        self.fint = zeros(shape=(nDof))
        self.mass = zeros(shape=(nDof, nDof))
        self.lumped = zeros(shape=(nDof))
        self.diss = 0.0

        self.outlabel = []

    def __str__(self):
        return self.state
