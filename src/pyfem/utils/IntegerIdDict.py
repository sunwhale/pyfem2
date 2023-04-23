from typing import List, Union


class IntegerIdDict(dict):
    """
    A list-like dictionary for storing and retrieving items with unique integer IDs.
    dict[key] <-> list(self.keys())
    id -> key
    indices -> list(self.keys())
    """

    def add_item_by_id(self, id_: int, item: object) -> None:
        """
        Add an item to the list-like dict with the specified ID.

        :param id_: the ID for the item
        :param item: the object to add
        :raises ValueError: if the ID already exists in the list
        """

        if id_ in self:
            raise ValueError(f"{type(self).__name__} already contains ID {id_}")

        self[id_] = item

    def get_items_by_ids(self, ids: Union[int, List[int]]) -> Union[object, List[object]]:
        """
        Get one or more items from the list-like dict by ID(s).

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
            raise TypeError("Argument to get_items_by_ids() must be int or list of ints")

    def get_indices_by_ids(self, ids: Union[int, List[int]] = None) -> Union[int, List[int]]:
        """
        Consider the keys of dict as a list: list(self.keys()).
        Get the indices of the list by the ids of the dict.

        :param ids: a single ID or a list of IDs (optional, defaults to all IDs in the list)
        :raises TypeError: if the argument is not an int, list of ints, or None
        :raises ValueError: if an ID is not found
        :return: either a single index or a list of indices
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
            raise TypeError("Argument to get_indices_by_ids() must be int, list of ints, or None")

    def get_id_by_index(self, index: int) -> int:
        """
        Get the ID of the item at the specified index position in the list(self.keys()).

        :param index: the index position in the list(self.keys())
        :raises IndexError: if the index is out of range
        :return: the ID of the item
        """

        indices = list(self.keys())

        try:
            return indices[index]
        except IndexError:
            raise IndexError("Index out of range")
