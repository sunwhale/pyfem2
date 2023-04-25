from numpy import zeros

from pyfem.utils.logger import get_logger
from pyfem.fem.NodeSet import NodeSet

logger = get_logger()


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
    def __init__(self, dictionary: dict = None):
        self.__dict__.update(dictionary or {})

    def __str__(self):
        props_list = []
        for key, value in self.__dict__.items():
            if key.startswith('_'):
                continue
            props_list.append(f'{key}: {value}')
        return '\n'.join(props_list)

    def __iter__(self):
        return iter(self.__dict__.items())

    def store(self, key: str, val: object) -> None:
        """
        store 方法用于动态添加属性和值。如果属性名中包含点号 .，则表示这是一个嵌套属性名，需要在类的 __dict__ 属性中按照层级结构创建一个字典对象，并在最终嵌套层级上设置属性值。
        如果属性名不包含点号，则直接在类的 __dict__ 属性中添加属性和值。
        注意，在 store 方法中使用了字典的 setdefault 方法，当属性名不存在时，会创建一个空的字典作为属性值。
        """
        if "." in key:
            keys = key.split(".")
            obj = self
            for k in keys[:-1]:
                obj = obj.__dict__.setdefault(k, {})
            obj.__dict__[keys[-1]] = clean_variable(val)
        else:
            self.__dict__[key] = clean_variable(val)


class GlobalData(Properties):
    def __init__(self, nodes: NodeSet, elements: "ElementSet", dofs):
        Properties.__init__(self, {'nodes': nodes, 'elements': elements, 'dofs': dofs})

        number_of_dofs = dofs.get_number_of_dofs()

        self.state = zeros(number_of_dofs)
        self.dstate = zeros(number_of_dofs)
        self.fint = zeros(number_of_dofs)
        self.fhat = zeros(number_of_dofs)

        self.velo = zeros(number_of_dofs)
        self.acce = zeros(number_of_dofs)

        self.solver_status = elements.solver_status

        self.outputNames = []

    def read_from_file(self, file_name):

        logger.info("Reading external forces ......")

        fin = open(file_name)

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

                            self.fhat[self.dofs.get_dof_ids_by_type(node_id, dof_type)] = eval(b[1])

    def printNodes(self, file_name=None, inodes=None):

        if file_name is None:
            f = None
        else:
            f = open(file_name, "w")

        if inodes is None:
            inodes = list(self.nodes.keys())

        print('   Node | ', file=f, end=' ')

        for dof_type in self.dofs.dof_types:
            print("  %-10s" % dof_type, file=f, end=' ')

        if hasattr(self, 'fint'):
            for dof_type in self.dofs.dof_types:
                print(" fint-%-6s" % dof_type, file=f, end=' ')

        for name in self.outputNames:
            print(" %-11s" % name, file=f, end=' ')

        print(" ", file=f)
        print(('-' * 100), file=f)

        for node_id in inodes:
            print('  %4i  | ' % node_id, file=f, end=' ')
            for dof_type in self.dofs.dof_types:
                print(' %10.3e ' % self.state[self.dofs.get_dof_ids_by_type(node_id, dof_type)], file=f, end=' ')
            for dof_type in self.dofs.dof_types:
                print(' %10.3e ' % self.fint[self.dofs.get_dof_ids_by_type(node_id, dof_type)], file=f, end=' ')

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
        self.dstate = elDstate
        self.stiff = zeros(shape=(nDof, nDof))
        self.fint = zeros(shape=(nDof))
        self.mass = zeros(shape=(nDof, nDof))
        self.lumped = zeros(shape=(nDof))
        self.diss = 0.0

        self.outlabel = []

    def __str__(self):
        return self.state
