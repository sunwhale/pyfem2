import pickle

from pyfem.utils.BaseModule import BaseModule


class DataDump(BaseModule):

    def __init__(self, props, globdat):

        self.prefix = globdat.prefix
        self.extension = ".dump"
        self.lastOnly = False

        BaseModule.__init__(self, props)

        if not hasattr(props, "interval"):
            self.interval = 1

        if self.lastOnly:
            self.interval = 1

    # ------------------------------------------------------------------------------
    #
    # ------------------------------------------------------------------------------

    def run(self, props, globdat):

        cycle = globdat.solverStatus.cycle

        if cycle % self.interval == 0:
            data = {}
            data["props"] = props
            data["globdat"] = globdat

            if self.lastOnly:
                name = str(self.prefix + self.extension)
            else:
                name = str(self.prefix + "_" + str(cycle) + self.extension)

            pickle.dump(data, open(name, "wb"))
