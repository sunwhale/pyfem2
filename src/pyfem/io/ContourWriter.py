from pyfem.utils.BaseModule import BaseModule
from pyfem.utils.logger import get_logger

logger = get_logger()


class ContourWriter(BaseModule):

    def __init__(self, props, globdat):

        self.prefix = globdat.prefix
        self.interval = 1
        self.k = 0

        BaseModule.__init__(self, props)

        self.columndata = []



    def run(self, props, globdat):

        if not globdat.SolverStatus.cycle % self.interval == 0:
            return

        logger.info("Writing contour file ......\n")

        crd = globdat.nodes.get_node_coords(self.nodes[0])
        outfile = open(self.prefix + '-contour-' + str(self.k) + '.out', 'w')

        outfile.write('#Node  %-10s %-10s' % ('x-coor', 'y-coor'))

        if len(crd) == 3:
            outfile.write('%-10s ' % 'z-coor')

        for dof_type in globdat.dofs.dof_types:
            outfile.write('%-10s ' % dof_type)

        for name in globdat.outputNames:
            outfile.write('%-10s ' % name)

        outfile.write('\n')

        for iNod in self.nodes:
            crd = globdat.nodes.get_node_coords(iNod)
            outfile.write('%4i %10.3e %10.3e' % (iNod, crd[0], crd[1]))

            if len(crd) == 3:
                outfile.write(' %10.3e' % crd[2])

            for dof_type in globdat.dofs.dof_types:
                outfile.write(' %10.3e' % (globdat.state[globdat.dofs.get_dof_ids_by_type(iNod, dof_type)]))

            for name in globdat.outputNames:
                stress = globdat.getData(name, list(range(len(globdat.nodes))))
                outfile.write(' %10.3e' % stress[iNod])

            outfile.write('\n')

        outfile.close()

        self.k = self.k + 1
