"""Chemical reaction edit 'attypes', 'bond', 'coord', and 'free_atom' data."""
from typing import List, Any
import numpy as np

from molecule_synthesizer.models import file_data


class BondError(Exception):
    pass


class LongBondError(Exception):
    pass


class Fragment(object):

    def __init__(self, name:str = None) -> None:
        self.name = name
        self.attypes = ''
        self.long_bond_idx = []
        self.bond = []
        self.coord = []
        self.free_atom = []

        self.attypes_obj = file_data.AttypeData(self.name)
        self.bond_obj = file_data.BondData(self.name)
        self.coord_obj = file_data.CoordData(self.name)
        self.xyz_obj = file_data.XYZFile(self.name)
        self.pdb_obj = file_data.PDBFile(self.name)

    def load_data(self) -> None:
        self.attypes = self.attypes_obj.load_data()
        self.long_bond_idx, self.bond = self.bond_obj.load_data()
        self.coord = self.coord_obj.load_data()

    def remove_hydrogen(self) -> None:
        """Remove Hydrogen which is seemed to be added.
        Process:
            1st: Remove a bond with hydrogen which is long and seemed to be added from bond data.
            2nd: Remove a hydrogen data from attypes data.
            3rd: Arrange bond index of bond data.
            4th: Remove a coord data of removed hydrogen from coord data.
            5th: Arrange a bond index of free-atom data.

        :return: new_fragment_data(dict): Edited fragment data that shows attypes, bond, long bond, free atom, and coord.
        """

        remover = RemoveHydrogen(self)
        remover.cut_bond()
        remover.remove_hydrogen_from_attype()
        remover.arrange_bond_idx()
        remover.remove_hydrogen_from_coord()
        remover.arrange_free_atom_idx()


    # def save_data(self):
    #     contents_creator = view.FileContentsCreator(self)
    #     attypes_contents = contents_creator.get_attypes_contents()
    #     bond_contents = contents_creator.get_bond_contents()
    #     coord_contents = contents_creator.get_coord_contents()
    #     xyz_contents = contents_creator.get_xyz_file_contents()
    #     pdb_contents = contents_creator.get_pdb_file_contents()
    #
    #     self.attypes_obj.save_data(attypes_contents)
    #     self.bond_obj.save_data(bond_contents)
    #     self.coord_obj.save_data(coord_contents)
    #     self.xyz_obj.save_data(xyz_contents)
    #     self.pdb_obj.save_data(pdb_contents)


class RemoveHydrogen:
    """Edit attypes, bond, index of bond, coord, index of free_atom."""
    def __init__(self, fragment: Fragment) -> None:
        self.fragment = fragment
        self.bond_to_cut = self.find_bond_to_cut()
        self.hydrogen_to_remove = self.find_hydrogen_to_remove()

    def find_bond_to_cut(self) -> Any:
        """Return one bond to cut.

        Return:
            bond_to_cut(list): One bond you can cut like, [1, 12],
            which means bond between atom(index=1) and atom(index=12) should be cut.
        TODO 切断できる箇所が何パターンかあったとき、どこを切断するのか確認する。
        """

        bonds = self.fragment.bond
        if not self.fragment.long_bond_idx:
            raise LongBondError('This fragment has no long bond.')
        bond_to_cut_idx = self.fragment.long_bond_idx.pop(0)

        return bond_to_cut_idx, bonds[bond_to_cut_idx]

    def find_hydrogen_to_remove(self) -> int:
        attype_data = self.fragment.attypes
        bond_idx, bond_to_cut = self.bond_to_cut

        atom1_idx = bond_to_cut[0]
        atom2_idx = bond_to_cut[1]

        if attype_data[atom1_idx] == 'H':
            return atom1_idx
        elif attype_data[atom2_idx] == 'H':
            return atom2_idx
        else:
            raise BondError("\"Remove Hydrogen\" class can only cut bond include hydrogen.")

    def cut_bond(self) -> None:
        """Choice one bond to cut and remove hydrogen"""
        bond_idx, bond_to_cut = self.bond_to_cut
        del self.fragment.bond[bond_idx]
        self.fragment.free_atom += bond_to_cut
        self.arrange_long_bond_idx()

    def arrange_long_bond_idx(self) -> None:
        long_bond_idx_list = self.fragment.long_bond_idx
        bond_idx_to_cut = self.bond_to_cut[0]
        for i, long_bond_idx in enumerate(long_bond_idx_list):
            if long_bond_idx > bond_idx_to_cut:
                long_bond_idx_list[i] -= 1
        self.fragment.long_bond_idx = long_bond_idx_list

    def remove_hydrogen(self) -> None:

        self.remove_hydrogen_from_attype()
        self.arrange_bond_idx()

        self.remove_hydrogen_from_coord()
        self.arrange_free_atom_idx()

    def remove_hydrogen_from_attype(self) -> None:
        atom_idx = self.hydrogen_to_remove
        del self.fragment.attypes[atom_idx]

    def arrange_bond_idx(self) -> None:
        """Organize index of bond data"""
        bonds = self.fragment.bond
        atom_idx = self.hydrogen_to_remove
        for i, bond in enumerate(self.fragment.bond):
            for j, atom_num in enumerate(bond):
                if atom_num > atom_idx:
                    bonds[i][j] -= 1
        self.fragment.bond = bonds

    def remove_hydrogen_from_coord(self) -> None:
        atom_idx = self.hydrogen_to_remove
        del self.fragment.coord[atom_idx]

    def arrange_free_atom_idx(self) -> None:
        atom_idx = self.hydrogen_to_remove
        free_atom = self.fragment.free_atom
        for i, atom in enumerate(free_atom):
            if atom == atom_idx:
                del free_atom[i]
            if atom > atom_idx:
                free_atom[i] -= 1
        self.fragment.free_atom = free_atom


class Synthesis:
    def __init__(self, fragment1: Fragment, fragment2: Fragment) -> None:
        self.fragment1 = fragment1
        self.fragment2 = fragment2

        new_mol_name = f'{fragment1.name}_{fragment2.name}'
        self.new_molecule = Fragment(new_mol_name)

    def synthesize_attypes(self):
        self.new_molecule.attypes = self.fragment1.attypes + self.fragment2.attypes[1:]

    def synthesize_bond(self):
        fragment1_element_num = len(self.fragment1.attypes) - 1
        updated_list = []
        for bond in self.fragment2.bond:
            updated_list.append([i + fragment1_element_num for i in bond])
        self.new_molecule.bond = self.fragment1.bond + updated_list

    def synthesize_free_atom(self):
        fragment1_element_num = len(self.fragment1.attypes) - 1
        adding_fragment_free_atom = [i + fragment1_element_num for i in self.fragment2.free_atom]
        self.new_molecule.free_atom = self.fragment1.free_atom + adding_fragment_free_atom


    def synthesize_long_bond(self):
        fragment1_bond_num = len(self.fragment1.bond)
        adding_fragment_long_bond_idx = [i + fragment1_bond_num for i in self.fragment2.long_bond_idx]
        self.new_molecule.long_bond_idx = self.fragment1.long_bond_idx + adding_fragment_long_bond_idx


    @staticmethod
    def get_fragment_x_min_max(coord_data: List[List[float]]) -> List[float]:
        """
        Return the maximum and minimum of each of the x, y, z of the box in which the molecule just fits.
        :return(list): [[x_min, x_max], [y_min, y_max], [z_min, z_max]]
        """
        coord_data = np.array(coord_data[1:])
        x_min, y_min, z_min = np.min(coord_data, axis=0)
        x_max, y_max, z_max = np.max(coord_data, axis=0)
        return [x_min, x_max]

    def get_x_diff(self) -> float:
        f1_max_x = self.get_fragment_x_min_max(self.fragment1.coord)[1]
        f2_min_x = self.get_fragment_x_min_max(self.fragment2.coord)[0]
        if f1_max_x < f2_min_x:
            return 0
        diff = f1_max_x - f2_min_x
        return diff

    def synthesize_coord(self):
        diff = self.get_x_diff()
        coord_f1 = self.fragment1.coord

        coord_f2 = np.array(self.fragment2.coord[1:])
        coord_f2 += np.array([diff+1, 0, 0])
        coord_f2 = list(coord_f2)
        self.new_molecule.coord = coord_f1 + coord_f2

    def add_bond(self):
        free_atom_1 = self.new_molecule.free_atom.pop(0)
        free_atom_2 = self.new_molecule.free_atom.pop(0)
        new_bond = [free_atom_1, free_atom_2]
        self.new_molecule.bond.append(new_bond)
