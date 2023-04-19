from numpy import zeros

from pyfem.fem.Assembly import assembleInternalForce, assembleExternalForce, assembleTangentStiffness, commit
from pyfem.utils.BaseModule import BaseModule
from pyfem.utils.logger import getLogger

logger = getLogger()




class LinearSolver(BaseModule):

    def __init__(self, props, globdat):
        BaseModule.__init__(self, props)

        self.fext = zeros(len(globdat.dofs))

        logger.info("Starting linear solver .......")



    def run(self, props, globdat):
        globdat.SolverStatus.increaseStep()

        K, fint = assembleTangentStiffness(props, globdat)
        fext = assembleExternalForce(props, globdat)

        state0 = globdat.state

        globdat.state = globdat.dofs.solve(K, fext)

        globdat.Dstate = globdat.state - state0

        globdat.fint = assembleInternalForce(props, globdat)

        commit(props, globdat)

        globdat.elements.commitHistory()

        globdat.active = False
