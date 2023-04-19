class BaseModule:

    def __init__(self, props):

        if hasattr(props, 'currentModule') and hasattr(props, props.currentModule):
            currentModule = props.currentModule
        else:
            currentModule = self.__class__.__name__

            if currentModule.endswith("olver") == "olver":
                currentModule = "solver"

        if hasattr(props, currentModule):
            self.myProps = getattr(props, currentModule)

            for name, val in self.myProps:
                setattr(self, name, val)
