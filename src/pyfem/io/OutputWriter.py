from pyfem.utils.BaseModule import BaseModule
from pyfem.utils.logger import get_logger

logger = get_logger()


class OutputWriter(BaseModule):

    def __init__(self, props, globdat):

        self.prefix = globdat.prefix + "_glob"
        self.extension = ".out"
        self.onScreen = False

        BaseModule.__init__(self, props)

        if not hasattr(self, "filename"):
            self.filename = self.prefix + self.extension



    def run(self, props, globdat):

        logger.info("Writing output file ..........")

        if self.onScreen:
            globdat.printNodes()

        globdat.printNodes(self.filename)
