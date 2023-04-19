class itemList(dict):

    def add(self, ID, item):

        if ID in self:
            raise RuntimeError('ID ' + str(ID) + ' already exists in ' + type(self).__name__)

        self[ID] = item

    def get(self, IDs):

        if isinstance(IDs, int):
            return self[IDs]
        elif isinstance(IDs, list):
            return [self[ID] for ID in IDs]

        raise RuntimeError('illegal argument for itemList.get')

    def getIndices(self, IDs=-1):

        if IDs == -1:
            return list(self.keys())
        elif isinstance(IDs, int):
            return list(self.keys()).index(IDs)
        elif isinstance(IDs, list):
            return [list(self.keys()).index(ID) for ID in IDs]

        raise RuntimeError('illegal argument for itemList.getIndices')

    def findID(self, index):

        return list(self.keys())[index]
