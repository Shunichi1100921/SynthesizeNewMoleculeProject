import collections
import os
from typing import List

from molecule_synthesizer.models.fragment import Fragment


class FileContentsCreator(object):
    """Create content from dictionary data for storage in a file."""
    def __init__(self, fragment: Fragment) -> None:
        self.fragment = fragment

    def get_attypes_contents(self) -> List[str]:
        """
        :param attype_data: List of attypes, the beginning of the list is blank.
                    ex) ['', 'C', 'O', 'H', 'H', 'O', 'C', 'H', 'H']
        :return: Insert a newline code (os.linesep) at the end of each element of the list.
        """
        attype_data = self.fragment.attypes
        contents = [f'{attype}{os.linesep}' for attype in attype_data[1:]]
        return contents

    def get_bond_contents(self) -> List[str]:
        """
        :param bond_data: List of bond list of atom num.
                    ex) [[1, 2], [1, 3], [4, 6], [5, 6], [5, 7], [6, 8], [1, 6]]
        :return: The atoms of each bond's are joined by a space, and a newline code is inserted at the end of each bond.
                    ex) ['1 2/n', '1 3/n', '4 6/n', '5 6/n', '5 7/n', '6 8/n', '1 6/n']
        """

        bond_data = self.fragment.bond
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

        coord_data = self.fragment.coord
        coord_str = self.list_of_list_of_num_to_list_of_str(coord_data[1:])
        contents = [f'{coord}{os.linesep}' for coord in coord_str]
        return contents

    def get_xyz_file_contents(self) -> List[str]:
        name = self.fragment.name
        attype_data = self.fragment.attypes[1:]
        coord_data = self.fragment.coord

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

    def get_pdb_file_contents(self) -> List[str]:
        hetatm = self.get_pdb_hetatm_contents()
        conect = self.get_pdb_conect_contents()
        contents = hetatm + conect
        contents.append('END\n')
        contents.insert(0, 'COMPND\nCOMPND\n')
        return contents

    def get_pdb_hetatm_contents(self) -> List[str]:
        attype = self.fragment.attypes
        coord = self.fragment.coord
        num_for_each_atom = collections.defaultdict(int)
        contents = []
        for i, atom in enumerate(attype):
            if not atom:
                continue
            num_for_each_atom[atom] += 1
            atom_index = num_for_each_atom[atom]
            coord_of_this_atom = coord[i]
            coord_of_this_atom = [round(i, 3) for i in coord_of_this_atom]

            contents_coord_part = f'{" "*(8-len(str(coord_of_this_atom[0])))}{coord_of_this_atom[0]}' \
                         f'{" "*(8-len(str(coord_of_this_atom[1])))}{coord_of_this_atom[1]}' \
                         f'{" "*(8-len(str(coord_of_this_atom[2])))}{coord_of_this_atom[2]}'
            contents_line = f'HETATM{" "*(5-len(str(i)))}{i} {atom}{atom_index}{" "*(17-len(str(atom_index)))}' \
                            f'{contents_coord_part}  1.00  0.00{" "*11}{atom}\n'
            contents.append(contents_line)
        return contents

    def get_pdb_conect_contents(self) -> List[str]:
        attype = self.fragment.attypes[1:]
        bonds = self.fragment.bond
        element_num = len(attype)
        contents = []
        for target_atom in range(1, element_num):
            data = self.create_pdb_conect_data(target_atom, bonds)
            contents_line = f'CONECT{data}\n'
            contents.append(contents_line)
        return contents

    def create_pdb_conect_data(self, target_atom: int, bonds: List[List[int]]) -> str:
        atom_list_connecting_target_include_target = []
        for bond in bonds:
            if target_atom in bond:
                atom_list_connecting_target_include_target += bond
        atom_list_connecting_target = [i for i in atom_list_connecting_target_include_target if not i == target_atom]
        contents_of_conect_part = f'{" "*(5-len(str(target_atom)))}{target_atom}'
        for atom in atom_list_connecting_target:
            contents_of_conect_part += f'{" "*(5-len(str(atom)))}{atom}'
        return contents_of_conect_part
