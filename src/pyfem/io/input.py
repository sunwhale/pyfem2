import argparse
import os.path
import pickle
import sys

from pyfem.fem.Contact import Contact
from pyfem.fem.DofSpace import DofSpace
from pyfem.fem.ElementSet import ElementSet
from pyfem.fem.NodeSet import NodeSet
from pyfem.utils.data_structures import Properties, GlobalData
from pyfem.utils.logger import set_logger
from pyfem.utils.parser import file_parser


def input_reader():
    inp_file_name, out_file_name, parameters = get_arguments()

    props = Properties()

    if out_file_name is not None:
        with open(out_file_name, 'rb') as f:
            data = pickle.load(f)
            props = data["props"]

    if inp_file_name is not None:
        if inp_file_name[-4:] == '.pro':
            props = file_parser(inp_file_name)
        else:
            props = file_parser(inp_file_name + '.pro')

    if parameters is not None:
        for p in parameters:
            x = p.split("=")
            props.store(x[0], x[1])

    if out_file_name is not None:
        return props, data["globdat"]

    input_file_name = props.input

    logger = set_logger(props)

    nodes = NodeSet()
    nodes.read_from_file(input_file_name)
    logger.info(nodes)

    elems = ElementSet(nodes, props)
    elems.read_from_file(input_file_name)
    logger.info(elems)

    dofs = DofSpace(elems)
    dofs.read_from_file(input_file_name)

    globdat = GlobalData(nodes, elems, dofs)

    globdat.read_from_file(input_file_name)

    globdat.active = True
    globdat.prefix = os.path.splitext(inp_file_name)[0]

    globdat.contact = Contact(props)

    globdat.print_nodes(node_ids=[0])

    return props, globdat


def get_arguments() -> tuple[str, str, list[str]]:

    # 创建一个 argparse 解析器对象
    parser = argparse.ArgumentParser(add_help=False)

    # 添加程序输入文件选项
    parser.add_argument('-i', metavar='input', type=str,
                        help='Identify the input file name.')

    # 添加程序输出文件选项
    parser.add_argument('-o', metavar='output', type=str,
                        help='Identify the output file name.')

    # 添加参数选项
    parser.add_argument('-p', metavar='parameter', type=str, action='append',
                        help='Parameters to pass to the program.')

    # 添加帮助选项
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                        help='Show this help message and exit.')

    # 添加版本选项
    parser.add_argument('-v', '--version', action='version', help='Show program\'s version number and exit.',
                        version='%(prog)s 0.1')

    # 解析命令行参数
    args = parser.parse_args()

    # 如果未指定程序输入文件，则打印帮助并退出
    if not args.i:
        print('error: the input file is required.')
        parser.print_help()
        sys.exit()

    return args.i, args.o, args.p
