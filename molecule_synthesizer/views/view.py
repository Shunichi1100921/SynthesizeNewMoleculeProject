"""Create file contents."""
from __future__ import annotations

import collections
import os

from molecule_synthesizer.models.fragment import Fragment


class FileContentsCreator(object):
    """Create content from dictionary data for storage in a file."""
    def __init__(self, fragment: Fragment) -> None:
        self.fragment = fragment

    def get_attypes_contents(self) -> list[str]:
        """Remove elements for index alignment and insert line a newline code at the end of each element of the list.

        Returns:
            contents (list): Contents list for 'attypes.dat' file.

        """
        attype_data = self.fragment.attypes
        contents = [f'{attype}{os.linesep}' for attype in attype_data[1:]]
        return contents

    def get_bond_contents(self) -> list[str]:
        """Change type int to type strings, separate the index of the atoms in each bond with a space,
        and insert a newline code at the end of each element of the list.

        Return:
             contents (list): Contents list for 'bond.dat' file.

        """

        bond_data = self.fragment.bond
        bond_str = self.list_of_list_of_num_to_list_of_str(bond_data)
        contents = [f'{bond}{os.linesep}' for bond in bond_str]
        return contents

    def get_coord_contents(self) -> list[str]:
        """Separate x-coordinate, y-coordinate, and z-coordinate with a space
        and insert a newline code at the end of each element of the list.

        Return:
            contents (list): Contents list for 'coord.dat' file.

        """

        coord_data = self.fragment.coord
        coord_str = self.list_of_list_of_num_to_list_of_str(coord_data[1:])
        contents = [f'{coord}{os.linesep}' for coord in coord_str]
        return contents

    def get_xyz_file_contents(self) -> list[str]:
        """The first two elements of the list are the number of atoms and the name of the molecule,
        followed by the atom type and the position of the atom, separated by spaces.

        Returns:
            contents (list): Contents list for '{molecule_name}.xyz' file.

        """
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

    def list_of_list_of_num_to_list_of_str(self, data: list[list[int]]) -> list[str]:
        """Utility function for change type of element of the list.

        ex) [[1, 2], [3, 4]] -> ['1 2', '3 4'], [[1.23, -2.13], [2,19, -12,4]] -> ['1.23 -2.13', '2,19 -12,4']

        Args:
            data: List of list of int/float.
        Return:
             List of str.
        """
        data_str = [list(map(str, l)) for l in data]
        return [' '.join(element) for element in data_str]

    def get_pdb_file_contents(self) -> list[str]:
        """Get pdb contents by dividing into two parts: HETATM, which is the part that represents the information of
        each atom, and CONECT, which represents the information of bindings.

        Returns:
            content (list): Contents for '{molecule_name}.pdb' file.

        """
        hetatm = self.get_pdb_hetatm_contents()
        conect = self.get_pdb_conect_contents()
        contents = hetatm + conect
        contents.append('END\n')
        contents.insert(0, 'COMPND\nCOMPND\n')
        return contents

    def get_pdb_hetatm_contents(self) -> list[str]:
        """Get "HETATM" contents.  It contains atom type and position.

        Returns:
            contents (list): HETATM part of pdb file.

        """
        attype = self.fragment.attypes
        coord = self.fragment.coord
        num_for_each_atom = collections.defaultdict(int)
        contents = []
        for i, atom in enumerate(attype):
            if not atom:
                continue
            num_for_each_atom[atom] += 1
            atom_index = num_for_each_atom[atom]
            atom_atom_idx = f'{atom}{atom_index}'
            coord_of_this_atom = coord[i]
            coord_of_this_atom = [round(i, 3) for i in coord_of_this_atom]

            contents_line = f'HETATM{i:>5} {atom_atom_idx:<3}{coord_of_this_atom[0]:>23}' \
                            f'{coord_of_this_atom[1]:>8}{coord_of_this_atom[2]:>8}  1.00 0.00{" "*11}{atom:<2}\n'

            contents.append(contents_line)
        return contents

    def get_pdb_conect_contents(self) -> list[str]:
        """Get "CONECT" contents.  It shows the atom index which each atom bind to.

        Returns:
            contents (list): "CONECT" part of pdb file.

        """
        attype = self.fragment.attypes[1:]
        element_num = len(attype)
        contents = []
        for target_atom in range(1, element_num):
            data = self.create_pdb_conect_data(target_atom)
            contents_line = f'CONECT{data}\n'
            contents.append(contents_line)
        return contents

    def create_pdb_conect_data(self, target_atom: int) -> str:
        """Return a string containing the indexes of the atoms bounded to the target atom.

        Args:
            target_atom (int): The index of the target atom.

        Returns:
            contents_of_conect_part (str): A string the index of the bonded atoms written in every four letters.

        """
        bonds = self.fragment.bond
        atom_list_connecting_target_include_target = []
        for bond in bonds:
            if target_atom in bond:
                atom_list_connecting_target_include_target += bond
        atom_list_connecting_target = [i for i in atom_list_connecting_target_include_target if not i == target_atom]
        contents_of_conect_part = f'{target_atom:>5}'
        for atom in atom_list_connecting_target:
            contents_of_conect_part += f'{atom:>5}'
        return contents_of_conect_part

    def get_fragid_contents(self) -> list[str]:
        """Get fragid contents.  It shows that which fragment each fragment is derived from.

        Returns:
            contents (list): Each string containing atom index and fragment type.

        """
        contents = [f'{i+1:<4}{fragid}\n' for i, fragid in enumerate(self.fragment.fragid)]
        return contents
