from numpy import zeros

from pyfem.fem.Assembly import assembleInternalForce, assembleExternalForce, assembleTangentStiffness, commit
from pyfem.utils.BaseModule import BaseModule
from pyfem.utils.logger import get_logger

logger = get_logger()




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

        globdat.dstate = globdat.state - state0

        globdat.fint = assembleInternalForce(props, globdat)

        commit(props, globdat)

        globdat.elements.commit_history()

        globdat.active = False
