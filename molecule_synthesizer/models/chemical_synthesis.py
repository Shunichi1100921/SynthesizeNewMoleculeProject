"""Chemical reaction edit 'attypes', 'bond', 'coord', and 'free_atom' data."""
from typing import List
import numpy as np


class BondError(Exception):
    pass


class LongBondError(Exception):
    pass


class RemoveHydrogen:
    """Edit attypes, bond, index of bond, coord, index of free_atom."""
    def __init__(self, fragment_data):
        self.fragment_data = fragment_data
        self.name = fragment_data['name']
        self.attypes = fragment_data['attypes']
        self.long_bond_idx = fragment_data['long_bond']
        self.bond = fragment_data['bond']
        self.coord = fragment_data['coord']
        self.free_atom = fragment_data['free_atom']
        self.bond_to_cut = self.find_bond_to_cut()
        self.hydrogen_to_remove = self.find_hydrogen_to_remove()

    def find_bond_to_cut(self):
        """Return one bond to cut.

        Return:
            bond_to_cut(list): One bond you can cut like, [1, 12],
            which means bond between atom(index=1) and atom(index=12) should be cut.
        TODO 切断できる箇所が何パターンかあったとき、どこを切断するのか確認する。
        """

        bonds = self.bond
        if not self.long_bond_idx:
            raise LongBondError('This fragment has no long bond.')
        bond_to_cut_idx = self.long_bond_idx.pop(0)

        return bond_to_cut_idx, bonds[bond_to_cut_idx]

    def find_hydrogen_to_remove(self):
        attype_data = self.attypes
        bond_idx, bond_to_cut = self.bond_to_cut

        atom1_idx = bond_to_cut[0]
        atom2_idx = bond_to_cut[1]

        if attype_data[atom1_idx] == 'H':
            return atom1_idx
        elif attype_data[atom2_idx] == 'H':
            return atom2_idx
        else:
            raise BondError("\"Remove Hydrogen\" class can only cut bond include hydrogen.")

    def cut_bond(self):
        """Choice one bond to cut and remove hydrogen"""
        bond_idx, bond_to_cut = self.bond_to_cut
        del self.bond[bond_idx]
        self.free_atom += bond_to_cut
        self.arrange_long_bond_idx()

    def arrange_long_bond_idx(self):
        long_bond_idx_list = self.long_bond_idx
        bond_idx_to_cut = self.bond_to_cut[0]
        for i, long_bond_idx in enumerate(long_bond_idx_list):
            if long_bond_idx > bond_idx_to_cut:
                long_bond_idx_list[i] -= 1
        self.long_bond_idx = long_bond_idx_list

    def remove_hydrogen(self):

        self.remove_hydrogen_from_attype()
        self.arrange_bond_idx()

        self.remove_hydrogen_from_coord()
        self.arrange_free_atom_idx()

    def remove_hydrogen_from_attype(self):
        atom_idx = self.hydrogen_to_remove
        del self.attypes[atom_idx]

    def arrange_bond_idx(self):
        """Organize index of bond data"""
        bonds = self.bond
        atom_idx = self.hydrogen_to_remove
        for i, bond in enumerate(self.bond):
            for j, atom_num in enumerate(bond):
                if atom_num > atom_idx:
                    bonds[i][j] -= 1
        self.bond = bonds

    def remove_hydrogen_from_coord(self):
        atom_idx = self.hydrogen_to_remove
        del self.coord[atom_idx]

    def arrange_free_atom_idx(self):
        atom_idx = self.hydrogen_to_remove
        free_atom = self.free_atom
        for i, atom in enumerate(free_atom):
            if atom == atom_idx:
                del free_atom[i]
            if atom > atom_idx:
                free_atom[i] -= 1
        self.free_atom = free_atom


class Synthesis:
    def __init__(self, fragment1_data, fragment2_data):
        self.fragment1_data = fragment1_data
        self.fragment2_data = fragment2_data

        self.name_f1 = fragment1_data['name']
        self.attypes_f1 = fragment1_data['attypes']
        self.long_bond_idx_f1 = fragment1_data['long_bond']
        self.bond_f1 = fragment1_data['bond']
        self.coord_f1 = fragment1_data['coord']
        self.free_atom_f1 = fragment1_data['free_atom']

        self.name_f2 = fragment2_data['name']
        self.attypes_f2 = fragment2_data['attypes']
        self.long_bond_idx_f2 = fragment2_data['long_bond']
        self.bond_f2 = fragment2_data['bond']
        self.coord_f2 = fragment2_data['coord']
        self.free_atom_f2 = fragment2_data['free_atom']

        self.new_mol_data = {}
        self.name_new_mol = ''
        self.attypes_new_mol = []
        self.long_bond_idx_new_mol = []
        self.bond_new_mol = []
        self.coord_new_mol = []
        self.free_atom_new_mol = []

    def synthesize_name(self):
        self.name_new_mol = f'{self.name_f1}_{self.name_f2}'

    def synthesize_attypes(self):
        self.attypes_new_mol = self.attypes_f1 + self.attypes_f2[1:]

    def synthesize_coord(self):
        self.coord_new_mol = self.coord_f1 + self.coord_f2[1:]

    def arrange_atom_index_of_bond(self):
        fragment1_element_num = len(self.attypes_f1) - 1
        updated_list = []
        for bond in self.bond_f2:
            updated_list.append([i + fragment1_element_num for i in bond])
        self.bond_f2 = updated_list

    def arrange_atom_index_of_free_atom(self):
        fragment1_element_num = len(self.attypes_f1) - 1
        self.free_atom_f2 = [i + fragment1_element_num for i in self.free_atom_f2]

    def synthesize_bond(self):
        self.bond_new_mol = self.bond_f1 + self.bond_f2

    def synthesize_free_atom(self):
        self.free_atom_new_mol = self.free_atom_f1 + self.free_atom_f2

    def arrange_bond_index_of_long_bond(self):
        fragment1_bond_num = len(self.bond_f1)
        self.long_bond_idx_f2 = [i + fragment1_bond_num for i in self.long_bond_idx_f2]

    def synthesize_long_bond(self):
        self.long_bond_idx_new_mol = self.long_bond_idx_f1 + self.long_bond_idx_f2

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
        coord_f1 = self.coord_f1
        coord_f2 = self.coord_f2
        f1_max_x = self.get_fragment_x_min_max(coord_f1)[1]
        f2_min_x = self.get_fragment_x_min_max(coord_f2)[0]
        if f1_max_x < f2_min_x:
            return -1
        diff = f1_max_x - f2_min_x
        return diff

    def arrange_coord_to_place_fragment(self):
        diff = self.get_x_diff()
        coord_f1 = self.coord_f1
        coord_f2 = np.array(self.coord_f2[1:])
        coord_f2 += np.array([diff+1, 0, 0])
        coord_f2 = coord_f2.tolist()
        self.coord_new_mol = coord_f1 + coord_f2

    def add_bond(self):
        free_atom_1 = self.free_atom_new_mol.pop(0)
        free_atom_2 = self.free_atom_new_mol.pop(0)
        new_bond = [free_atom_1, free_atom_2]
        self.bond_new_mol.append(new_bond)
