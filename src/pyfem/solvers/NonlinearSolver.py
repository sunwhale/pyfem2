import sys

from numpy import zeros

from pyfem.fem.Assembly import assembleExternalForce
from pyfem.fem.Assembly import assembleTangentStiffness
from pyfem.utils.BaseModule import BaseModule
from pyfem.utils.logger import get_logger

logger = get_logger()


class NonlinearSolver(BaseModule):

    def __init__(self, props, globdat):

        self.tol = 1.0e-3
        self.iterMax = 10

        self.maxCycle = sys.maxsize
        self.maxLam = 1.0e20
        self.dtime = 1.0
        self.loadFunc = "t"
        self.loadCases = []

        BaseModule.__init__(self, props)

        if self.maxLam > 1.0e19 and self.maxCycle == sys.maxsize:
            self.maxCycle = 5

        globdat.lam = 0.0
        globdat.solver_status.dtime = self.dtime

        self.loadfunc = eval("lambda t : " + str(self.loadFunc))

        if hasattr(self, "loadTable"):
            self.maxCycle = len(self.loadTable)
            loadTable = zeros(self.maxCycle + 1)
            loadTable[1:] = self.loadTable
            self.loadTable = loadTable

        logger.info("Starting nonlinear solver .........")

    # ------------------------------------------------------------------------------
    #
    # ------------------------------------------------------------------------------

    def run(self, props, globdat):

        stat = globdat.solver_status

        stat.increaseStep()

        dofCount = globdat.dofs.number_of_dofs

        a = globdat.state
        Da = globdat.dstate

        Da[:] = zeros(dofCount)
        fint = zeros(dofCount)

        logger.info("Nonlinear solver ............")
        logger.info("    =============================================")
        logger.info("    Load step %i" % globdat.solver_status.cycle)
        logger.info("    =============================================")
        logger.info('    Newton-Raphson   : L2-norm residual')

        self.setLoadAndConstraints(globdat)

        K, fint = assembleTangentStiffness(props, globdat)

        error = 1.

        self.setLoadAndConstraints(globdat)

        fext = assembleExternalForce(props, globdat)

        while error > self.tol:

            stat.iiter += 1

            da = globdat.dofs.solve(K, fext - fint)

            Da[:] += da[:]
            a[:] += da[:]

            K, fint = assembleTangentStiffness(props, globdat)

            # note that the code is different from the one presented in the book, which
            # is slightly shorter for the sake of clarity.
            # In the case of a prescribed displacement, the external force is zero
            # and hence its norm is zero. In that case, the norm of the residue is not
            # divided by the norm of the external force.

            norm = globdat.dofs.norm(fext)

            if norm < 1.0e-16:
                error = globdat.dofs.norm(fext - fint)
            else:
                error = globdat.dofs.norm(fext - fint) / norm

            logger.info('    Iteration %4i   : %6.4e' % (stat.iiter, error))

            globdat.dofs.set_constrain_factor(0.0)

            if stat.iiter == self.iterMax:
                raise RuntimeError('Newton-Raphson iterations did not converge!')

        # Converged

        globdat.elements.update_commit_history()

        Da[:] = zeros(globdat.dofs.number_of_dofs)

        globdat.fint = fint

        if stat.cycle == self.maxCycle or globdat.lam > self.maxLam:
            globdat.active = False

        # -------------------------------------------------------------------------------

    #
    # -------------------------------------------------------------------------------

    def setLoadAndConstraints(self, globdat):

        if hasattr(self, "loadTable"):
            cycle = globdat.solver_status.cycle

            globdat.lam = self.loadTable[cycle]
            globdat.dlam = self.loadTable[cycle] - self.loadTable[cycle - 1]

            globdat.dofs.set_constrain_factor(globdat.dlam)

            globdat.solver_status.lam = globdat.lam
        else:
            globdat.lam = self.loadfunc(globdat.solver_status.time)
            lam0 = self.loadfunc(globdat.solver_status.time - globdat.solver_status.dtime)

            globdat.dlam = globdat.lam - lam0
            globdat.dofs.set_constrain_factor(globdat.dlam)

            globdat.solver_status.lam = globdat.lam

            logger.debug('  ---- main load -------------------------')
            logger.debug('    loadFactor       : %4.2f' % globdat.lam)
            logger.debug('    incr. loadFactor : %4.2f' % globdat.dlam)

            for loadCase in self.loadCases:
                loadProps = getattr(self.myProps, loadCase)
                loadfunc = eval("lambda t : " + str(loadProps.loadFunc))
                lam = loadfunc(globdat.solver_status.time)
                lam0 = loadfunc(globdat.solver_status.time - globdat.solver_status.dtime)
                dlam = lam - lam0
                globdat.dofs.set_constrain_factor(dlam, loadProps.nodeTable)

                logger.debug('  ---- %s ---------------------' % loadCase)
                logger.debug('    loadFactor       : %4.2f' % lam)
                logger.debug('    incr. loadFactor : %4.2f' % dlam)
