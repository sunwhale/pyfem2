from numpy import ndarray
from pylab import plot, xlabel, ylabel, ion, gcf

from pyfem.utils.BaseModule import BaseModule
from pyfem.utils.data_structures import Properties


class GraphWriter(BaseModule):

    def __init__(self, props, globdat):

        self.prefix = globdat.prefix
        self.extension = ".out"
        self.onScreen = False

        BaseModule.__init__(self, props)

        if not hasattr(self, "filename"):
            self.filename = self.prefix + self.extension

        self.columndata = []

        for i, col in enumerate(self.columns):

            if hasattr(self, col):
                colProps = getattr(self, col)
            else:
                colProps = Properties()

            if not hasattr(colProps, "type"):
                colProps.type = col

            if not hasattr(colProps, "factor"):
                colProps.factor = 1.0

            if hasattr(colProps, "node"):
                if type(colProps.node) == str:
                    colProps.node = globdat.nodes.groups[colProps.node]

            self.columndata.append(colProps)

        if self.onScreen and hasattr(globdat, "onScreen"):
            self.onScreen = False
        else:
            globdat.onScreen = True

            self.fig = gcf()
            self.fig.show()
            self.fig.canvas.draw()

        self.outfile = open(self.filename, 'w')

        if self.onScreen:
            self.output = []

            xlabel(self.columns[0])
            ylabel(self.columns[1])

            ion()
        else:
            self.output = None

        self.run(props, globdat)

    def run(self, props, globdat):

        a = []

        for i, col in enumerate(self.columndata):

            if col.type in globdat.outputNames:
                data = globdat.getData(col.type, col.node)

            elif hasattr(globdat, col.type):
                b = getattr(globdat, col.type)
                if type(b) is ndarray:
                    if type(col.node) is list:
                        data = 0.0
                        for nod in col.node:
                            data += b[globdat.dofs.get_dof_ids_by_type(int(nod), col.dof)]
                    else:
                        data = b[globdat.dofs.get_dof_ids_by_type(col.node, col.dof)]
                else:
                    data = b
            elif col.type in globdat.outputNames:
                data = globdat.getData(col.type, col.node)
            elif hasattr(globdat.SolverStatus, col.type):
                data = getattr(globdat.SolverStatus, col.type)

            data = data * col.factor

            a.append(data)

            self.outfile.write(str(data) + ' ', )
            self.outfile.flush()

        self.outfile.write('\n')

        self.output.append(a)

        plot([x[0] for x in self.output], [x[1] for x in self.output], 'ro-')
        # plot( [x[0] for x in self.output], [x[2] for x in self.output], 'bo-' )    

        if self.onScreen:
            self.fig.canvas.draw()

        self.fig.savefig(self.prefix + '.png')
        '''
        if self.onScreen: 
          self.output.append( a )
    
          plot( [x[0] for x in self.output], [x[1] for x in self.output], 'ro-' )
          self.fig.canvas.draw()
        '''

        if not globdat.active:
            self.outfile.close()
