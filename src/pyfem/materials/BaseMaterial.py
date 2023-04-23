import copy

from numpy import zeros


class BaseMaterial:

    def __init__(self, props):

        self.numericalTangent = False
        self.storeOutputFlag = False

        for name, val in props:
            setattr(self, name, val)

        self.oldHistory = {}
        self.newHistory = {}

        self.outLabels = []
        self.solver_status = props.solver_status

    def setHistoryParameter(self, name, val):

        self.newHistory[name] = val
        return

    def getHistoryParameter(self, name):

        if type(self.oldHistory[name]) == float:
            return self.oldHistory[name]
        else:
            return self.oldHistory[name].copy()

    def commit_history(self):

        self.oldHistory = copy.deepcopy(self.newHistory)

    def setOutputLabels(self, labels):

        self.outLabels = labels
        self.outData = zeros(len(self.outLabels))
        return

    def storeOutputs(self, data):
        if self.storeOutputFlag:
            self.outData = data
        return
