class Solver:

    def __init__(self, props, globdat):
        solverProps = getattr(props, "solver")

        solverType = solverProps.type

        exec("from pyfem.solvers." + solverType + " import " + solverType)

        props.currentModule = "solver"

        self.solver = eval(solverType + "( props , globdat )")

    # -------------------------------------------------------------------------------
    #
    # -------------------------------------------------------------------------------

    def run(self, props, globdat):
        self.solver.run(props, globdat)
