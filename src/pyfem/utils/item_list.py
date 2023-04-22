from typing import List, Union


class ItemList(dict):
    """
    A dictionary-like object for storing and retrieving items with unique integer IDs.
    """

    def add(self, id_: int, item: object) -> None:
        """
        Add an item to the list with the specified ID.

        :param id_: the ID for the item
        :param item: the item to add
        :raises ValueError: if the ID already exists in the list
        """

        if id_ in self:
            raise ValueError(f"{type(self).__name__} already contains ID {id_}")

        self[id_] = item

    def get(self, ids: Union[int, List[int]]) -> Union[object, List[object]]:
        """
        Get one or more items from the list by ID(s).

        :param ids: a single ID or a list of IDs
        :raises TypeError: if the argument is not an int or list of ints
        :raises KeyError: if an ID is not found
        :return: either a single item or a list of items
        """

        if isinstance(ids, int):
            return self[ids]
        elif isinstance(ids, list):
            return [self[id_] for id_ in ids]
        else:
            raise TypeError("Argument to get() must be int or list of ints")

    def get_indices(self, ids: Union[int, List[int]] = None) -> Union[int, List[int]]:
        """
        Get the index position of one or more IDs in the list.

        :param ids: a single ID or a list of IDs (optional, defaults to all IDs in the list)
        :raises TypeError: if the argument is not an int, list of ints, or None
        :raises ValueError: if an ID is not found in the list
        :return: either a single index, a list of indices, or a list of all indices in the list
        """

        indices = list(self.keys())

        if ids is None:
            return indices
        elif isinstance(ids, int):
            try:
                return indices.index(ids)
            except ValueError:
                raise ValueError(f"ID {ids} not found in {type(self).__name__}")
        elif isinstance(ids, list):
            return [i for i, id_ in enumerate(indices) if id_ in ids]
        else:
            raise TypeError("Argument to get_indices() must be int, list of ints, or None")

    def find_id(self, index: int) -> int:
        """
        Get the ID of the item at the specified index position in the list.

        :param index: the index position of the item
        :raises IndexError: if the index is out of range
        :return: the ID of the item
        """

        try:
            return list(self.keys())[index]
        except IndexError:
            raise IndexError("Index out of range")
