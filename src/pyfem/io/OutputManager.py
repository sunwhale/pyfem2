class OutputManager:

    def __init__(self, props, globdat):

        self.outman = []

        outputModules = props.outputModules

        for name in outputModules:

            props.currentModule = name

            ioType = name

            if hasattr(props, name):
                moduleProps = getattr(props, name)
                if hasattr(moduleProps, "type"):
                    ioType = moduleProps.type

            exec("from pyfem.io." + ioType + " import " + ioType)

            self.outman.append(eval(ioType + "( props , globdat )"))

    def run(self, props, globdat):

        for i, output in enumerate(self.outman):
            output.run(props, globdat)
