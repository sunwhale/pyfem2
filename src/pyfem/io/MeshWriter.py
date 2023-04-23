from pyfem.utils.BaseModule import BaseModule
from pyfem.utils.logger import get_logger

logger = get_logger()




class MeshWriter(BaseModule):

    def __init__(self, props, globdat):

        self.prefix = globdat.prefix
        self.elementGroup = "All"
        self.k = 0
        self.interval = 1
        self.extraFields = []
        self.beam = False
        self.interface = False

        BaseModule.__init__(self, props)

        if type(self.extraFields) is str:
            self.extraFields = [self.extraFields]

    def run(self, props, globdat):

        if not globdat.SolverStatus.cycle % self.interval == 0:
            return

        logger.info("Writing mesh .................")

        dim = globdat.state.ndim

        if dim == 1:
            self.writeCycle(globdat.state, props, globdat)
        elif dim == 2:
            for state in globdat.state.transpose():
                self.writeCycle(state, props, globdat)

        self.writePvd()

    #
    #
    #

    def writeCycle(self, state, props, globdat):

        vtkfile = open(self.prefix + '-' + str(self.k) + '.vtu', 'w')

        vtkfile.write('<?xml version="1.0"?>\n')
        vtkfile.write(
            '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
        vtkfile.write('<UnstructuredGrid>\n')
        vtkfile.write('<Piece NumberOfPoints="' + str(len(globdat.nodes)) + '" NumberOfCells="')
        vtkfile.write(str(globdat.elements.element_group_count(self.elementGroup)) + '">\n')
        vtkfile.write('<PointData>\n')
        vtkfile.write('<DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii" >\n')

        dispDofs = ["u", "v", "w"]

        for node_id in list(globdat.nodes.keys()):
            for dispDof in dispDofs:
                if dispDof in globdat.dofs.dof_types:
                    vtkfile.write(str(state[globdat.dofs.getForType(node_id, dispDof)]) + ' ')
                else:
                    vtkfile.write(' 0.\n')

        vtkfile.write('</DataArray>\n')

        for field in self.extraFields:
            vtkfile.write('<DataArray type="Float64" Name="' + field + '" NumberOfComponents="1" format="ascii" >\n')

            for node_id in list(globdat.nodes.keys()):
                vtkfile.write(str(state[globdat.dofs.getForType(node_id, field)]) + ' ')

            vtkfile.write('</DataArray>\n')

        for name in globdat.outputNames:
            stress = globdat.getData(name, list(range(len(globdat.nodes))))

            vtkfile.write('<DataArray type="Float64" Name="' + name + '" NumberOfComponents="1" format="ascii" >\n')
            for i in range(len(globdat.nodes)):
                vtkfile.write(str(stress[i]) + " \n")

            vtkfile.write('</DataArray>\n')

        vtkfile.write('</PointData>\n')
        vtkfile.write('<CellData>\n')
        vtkfile.write('</CellData>\n')
        vtkfile.write('<Points>\n')
        vtkfile.write('<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n')

        for node_id in list(globdat.nodes.keys()):
            crd = globdat.nodes.get_node_coords(node_id)
            if len(crd) == 2:
                vtkfile.write(str(crd[0]) + ' ' + str(crd[1]) + " 0.0\n")
            else:
                vtkfile.write(str(crd[0]) + ' ' + str(crd[1]) + ' ' + str(crd[2]) + "\n")

        vtkfile.write('</DataArray>\n')
        vtkfile.write('</Points>\n')
        vtkfile.write('<Cells>\n')
        vtkfile.write('<DataArray type="Int64" Name="connectivity" format="ascii">\n')

        # --Store elements-----------------------------

        rank = globdat.nodes.rank

        for element in globdat.elements.iter_element_group(self.elementGroup):
            el_nodes = globdat.nodes.get_indices_by_ids(element.getNodes())

            if rank == 2:
                if len(el_nodes) == 2 and element.family == "BEAM":
                    vtkfile.write(str(el_nodes[0]) + ' ' + str(el_nodes[1]))
                elif len(el_nodes) == 3 and element.family == "BEAM":
                    vtkfile.write(str(el_nodes[0]) + ' ' + str(el_nodes[2]))
                elif len(el_nodes) == 3 or (len(el_nodes) == 4 and not self.interface):
                    for node in el_nodes:
                        vtkfile.write(str(node) + ' ')
                elif len(el_nodes) == 4 and self.interface:
                    vtkfile.write(str(el_nodes[0]) + ' ' + str(el_nodes[1]) + ' ' + str(el_nodes[3]) + ' ' + str(
                        el_nodes[2]) + ' ')
                elif len(el_nodes) == 6 or len(el_nodes) == 8:
                    for node in el_nodes[::2]:
                        vtkfile.write(str(node) + ' ')

            elif rank == 3:
                if len(el_nodes) <= 8:
                    for node in el_nodes:
                        vtkfile.write(str(node) + ' ')

            vtkfile.write('\n')

        vtkfile.write('</DataArray>\n')
        vtkfile.write('<DataArray type="Int64" Name="offsets" format="ascii">\n')

        nTot = 0

        for i, element in enumerate(globdat.elements.iter_element_group(self.elementGroup)):
            num_element_nodes = len(globdat.nodes.get_indices_by_ids(element.getNodes()))

            if rank == 2 and num_element_nodes == 8:
                num_element_nodes = 4
            elif num_element_nodes == 3 and self.beam:
                num_element_nodes = 2

            nTot += num_element_nodes
            vtkfile.write(str(nTot) + '\n')

        vtkfile.write('</DataArray>\n')
        vtkfile.write('<DataArray type="UInt8" Name="types" format="ascii">\n')

        for element in globdat.elements.iter_element_group(self.elementGroup):
            num_element_nodes = len(globdat.nodes.get_indices_by_ids(element.getNodes()))

            if rank == 2:
                if num_element_nodes < 4 and self.beam:
                    vtkfile.write('3\n')
                elif num_element_nodes == 3 or num_element_nodes == 6:
                    vtkfile.write('5\n')
                else:
                    vtkfile.write('9\n')
            else:
                if num_element_nodes == 8:
                    vtkfile.write('12\n')
                elif num_element_nodes == 6:
                    vtkfile.write('13\n')
                elif num_element_nodes == 4:
                    vtkfile.write('10\n')
                elif num_element_nodes == 5:
                    vtkfile.write('14\n')

        vtkfile.write('</DataArray>\n')
        vtkfile.write('</Cells>\n')
        vtkfile.write('</Piece>\n')
        vtkfile.write('</UnstructuredGrid>\n')
        vtkfile.write('</VTKFile>\n')

        self.k = self.k + 1

    
    #  writePvd
    

    def writePvd(self):

        f = open(self.prefix + '.pvd', 'w')

        f.write("<VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>\n")
        f.write("<Collection>\n")

        for i in range(self.k):
            f.write("<DataSet file='" + self.prefix + '-' + str(i) + ".vtu' groups='' part='0' timestep='" + str(
                i) + "'/>\n")

        f.write("</Collection>\n")
        f.write("</VTKFile>\n")

        f.close()
