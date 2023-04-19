import getopt
import os.path
import pickle

from pyfem.fem.Contact import Contact
from pyfem.fem.DofSpace import DofSpace
from pyfem.fem.ElementSet import ElementSet
from pyfem.fem.NodeSet import NodeSet
from pyfem.utils.dataStructures import GlobalData
from pyfem.utils.logger import setLogger
from pyfem.utils.parser import fileParser


def input_reader(argv):
    pName, dName, params = get_arguments(argv)
    return input_read(pName, dName, params)


def input_read(fname, dname=None, parameters=None):
    if dname is not None:
        with open(dname, 'rb') as f:
            data = pickle.load(f)
            props = data["props"]

    if fname is not None:
        if fname[-4:] == '.pro':
            props = fileParser(fname)
        else:
            props = fileParser(fname + '.pro')

    if parameters is not None:
        for p in parameters:
            x = p.split("=")
            props.store(x[0], x[1])

    if dname is not None:
        return props, data["globdat"]

    dataFileName = props.input

    logger = setLogger(props)

    nodes = NodeSet()
    nodes.readFromFile(dataFileName)
    logger.info(nodes)

    elems = ElementSet(nodes, props)
    elems.readFromFile(dataFileName)
    logger.info(elems)

    dofs = DofSpace(elems)
    dofs.readFromFile(dataFileName)

    globdat = GlobalData(nodes, elems, dofs)

    globdat.readFromFile(dataFileName)

    globdat.active = True
    globdat.prefix = os.path.splitext(fname)[0]

    globdat.contact = Contact(props)

    return props, globdat


def get_arguments(argv):
    slist = 'd:i:hvp:'
    llist = ['dump=', 'input=', 'help', 'version']

    options, remainder = getopt.getopt(argv[1:], slist, llist)

    pro_file_name = None
    dump_file_name = None
    parameters = []

    if len(options) == 0:
        pro_file_name = argv[1]
        options, remainder = getopt.getopt(argv[2:], slist, llist)

    for opt, arg in options:
        if opt in ('-i', '--input'):
            pro_file_name = arg
        elif opt in ('-d', '--dump'):
            dump_file_name = arg
        elif opt in ('-h', '--help'):
            print("Help")
        elif opt in ('-p', '--param'):
            parameters.append(arg)

    return pro_file_name, dump_file_name, parameters
