import os
from typing import List
from typing import Optional


class FileContentsCreator(object):
    """Create content from dictionary data for storage in a file."""
    def __init__(self, fragment_data):
        self.fragment_data = fragment_data
        self.name = fragment_data['name']
        self.attypes = fragment_data['attypes']
        self.coord = fragment_data['coord']
        self.bond = fragment_data['bond']

    def get_attypes_contents(self) -> List[str]:
        """
        :param attype_data: List of attypes, the beginning of the list is blank.
                    ex) ['', 'C', 'O', 'H', 'H', 'O', 'C', 'H', 'H']
        :return: Insert a newline code (os.linesep) at the end of each element of the list.
        """
        attype_data = self.attypes
        contents = [f'{attype}{os.linesep}' for attype in attype_data[1:]]
        return contents

    def get_bond_contents(self) -> List[str]:
        """
        :param bond_data: List of bond list of atom num.
                    ex) [[1, 2], [1, 3], [4, 6], [5, 6], [5, 7], [6, 8], [1, 6]]
        :return: The atoms of each bond's are joined by a space, and a newline code is inserted at the end of each bond.
                    ex) ['1 2/n', '1 3/n', '4 6/n', '5 6/n', '5 7/n', '6 8/n', '1 6/n']
        """

        bond_data = self.bond
        bond_str = self.list_of_list_of_num_to_list_of_str(bond_data)
        contents = [f'{bond}{os.linesep}' for bond in bond_str]
        return contents

    def get_coord_contents(self) -> List[str]:
        """
        :param coord_data: List of coord list of position on x, y, and z axes, the beginning of the list is blank.
                    ex) [[], [-0.4042, -1.0549, -0.7103], [0.5434, -1.1772, 0.2523], [0.2746, -1.4996, 1.4339]]

        :return: The atoms of each bond's are joined by a space, and a newline code is inserted at the end of each bond.
                    ex) ['-0.4042 -1.0549 -0.7103/n', '0.5434 -1.1772 0.2523/n', '0.2746 -1.4996 1.4339/n']
        """

        coord_data = self.coord
        coord_str = self.list_of_list_of_num_to_list_of_str(coord_data[1:])
        contents = [f'{coord}{os.linesep}' for coord in coord_str]
        return contents

    def get_xyz_file_contents(self) -> List[str]:
        name = self.name
        attype_data = self.attypes[1:]
        coord_data = self.coord

        coord_data = self.list_of_list_of_num_to_list_of_str(coord_data[1:])
        element_num = str(len(attype_data))

        data_list = [f'{attype} {coord}' for attype, coord in zip(attype_data, coord_data)]
        data_list.insert(0, element_num)
        data_list.insert(1, name)
        contents = [f'{data}{os.linesep}' for data in data_list]

        return contents

    def list_of_list_of_num_to_list_of_str(self, data) -> List[str]:
        """
        ex) [[1, 2], [3, 4]] -> ['1 2', '3 4'], [[1.23, -2.13], [2,19, -12,4]] -> ['1.23 -2.13', '2,19 -12,4']
        :param data: List of list of int/float.
        :return: List of str.
        """
        data_str = [list(map(str, l)) for l in data]
        return [' '.join(element) for element in data_str]
