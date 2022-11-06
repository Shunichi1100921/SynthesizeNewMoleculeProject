"""Chemical reaction edit 'attypes', 'bond', 'coord', and 'free_atom' data."""
from collections import defaultdict
from typing import List, Dict
import numpy as np

from molecule_synthesizer.models import file_data


class BondError(Exception):
    pass


class LongBondError(Exception):
    pass


class BindFragmentError(Exception):
    pass


class LossFragmentError(Exception):
    pass


class Fragment(object):

    def __init__(self, name: str = None, new_molecule: bool = False) -> None:
        self.name = name

        self.attypes_obj = file_data.AttypeData(self.name)
        self.bond_obj = file_data.BondData(self.name)
        self.coord_obj = file_data.CoordData(self.name)
        self.bind_fragment_obj = file_data.BindFragment(self.name)
        self.xyz_obj = file_data.XYZFile(self.name)
        self.pdb_obj = file_data.PDBFile(self.name)

        if new_molecule:
            self.attypes = []
            self.long_bond_idx = []
            self.bond = []
            self.coord = []
            self.bind_fragment = []
        else:
            self.load_data()

        self.free_atom = []

    def __len__(self) -> int:
        return len(self.attypes) - 1

    def load_data(self) -> None:
        self.attypes = self.attypes_obj.load_data()
        self.long_bond_idx, self.bond = self.bond_obj.load_data()
        self.coord = self.coord_obj.load_data()
        self.bind_fragment = self.bind_fragment_obj.load_data()

    def remove_hydrogen(self, bind_fragment: str) -> None:
        """Remove Hydrogen which is seemed to be added.
        Process:
            1st: Remove a bond with hydrogen which is long and seemed to be added from bond data.
            2nd: Remove a hydrogen data from attypes data.
            3rd: Arrange bond index of bond data.
            4th: Remove a coord data of removed hydrogen from coord data.
            5th: Arrange a bond index of free-atom data.

        """

        remover = RemoveHydrogen(self, bind_fragment)
        remover.remove_bond()


class RemoveHydrogen(object):

    def __init__(self, fragment: Fragment, bind_fragment: str) -> None:
        self.fragment = fragment
        self.bind_fragment = bind_fragment
        self.bond_idx_to_remove = self.find_bond_idx_to_remove()
        self.atom_idx_by_attype = self.find_atom_idx_by_attype()

    def find_bond_idx_to_remove(self) -> int:
        """Return one bond to cut.

        Return:
            bond idx to remove.
        """

        for bond_idx in self.fragment.long_bond_idx:
            bond = self.fragment.bond[bond_idx]
            for atom_idx in bond:
                if self.fragment.bind_fragment[atom_idx] == self.bind_fragment:
                    return bond_idx

        raise BindFragmentError("There is not appropriate fragment.")

    def find_atom_idx_by_attype(self) -> Dict[str, int]:
        """Look at the atoms of the bond to be cut and return hydrogen or bind_fragment.
        Return: {'hydrogen': atom_idx, 'bind_fragment': atom_idx}
            ex) {'hydrogen': 2, 'benzothiazole': 10}
        """
        attype_in_bond = {}
        for atom_idx in self.fragment.bond[self.bond_idx_to_remove]:
            if self.fragment.attypes[atom_idx] == 'H':
                attype_in_bond['hydrogen'] = atom_idx
            else:
                attype_in_bond[self.bind_fragment] = atom_idx

        return attype_in_bond

    def remove_bond(self) -> None:
        remove_atom_idx = self.atom_idx_by_attype['hydrogen']

        # Remove bond from self.fragment.bond.
        del self.fragment.bond[self.bond_idx_to_remove]

        # Add free atom.
        self.fragment.free_atom.append(self.atom_idx_by_attype[self.bind_fragment])

        # Remove hydrogen from self.fragment.attypes.
        del self.fragment.attypes[remove_atom_idx]

        # Arrange self.fragment.bond considering atom_idx.
        for i, bond in enumerate(self.fragment.bond):
            for j, atom_idx in enumerate(bond):
                if atom_idx > remove_atom_idx:
                    self.fragment.bond[i][j] -= 1

        # Remove hydrogen from self.fragment.coord.
        del self.fragment.coord[remove_atom_idx]

        # Arrange self.fragment.free_atom considering atom_idx.
        for i, atom_idx in enumerate(self.fragment.free_atom):
            if atom_idx > remove_atom_idx:
                self.fragment.free_atom[i] -= 1

        # Arrange self.fragment.long_bond_idx considering bond_idx.
        for i, bond_idx in enumerate(self.fragment.long_bond_idx):
            if bond_idx > self.bond_idx_to_remove:
                self.fragment.long_bond_idx[i] -= 1

        # Remove bond from self.fragment.long_bond_idx.
        self.fragment.long_bond_idx.remove(self.bond_idx_to_remove)

        # Remove hydrogen from self.fragment.bind_fragment.
        del self.fragment.bind_fragment[remove_atom_idx]


class Synthesis:
    """
    Synthesis Step:
        1st: synthesize attypes
        2nd: synthesize bond
            1st: synthesize bond
            2nd: arrange atom idx of bond
        3rd: synthesize free atom <- Don't need?
            1st: synthesize free atom
            2nd: arrange atom idx of bond
        4th: synthesize long_bond <- Don't need?
            1st: synthesize long_bond
            2nd: arrange bond idx of long_bond
        5th: synthesize coord
            1st: find appropriate position
            2nd: synthesize coord
        6th: add bond
            1st: synthesize atom idx of bond
            2nd: add bond
            3rd: remove free atom
        7th: synthesize bind_fragment
    """

    def __init__(self, fragments: Dict[str, str]) -> None:

        self.fragments = []

        # Make fragment object.
        if 'benzothiazole' not in fragments:
            raise LossFragmentError('Lack necessary fragment.')
        else:
            self.benzothiazole = Fragment(fragments['benzothiazole'])
            self.fragments.append(self.benzothiazole)

        if 'amide' not in fragments:
            raise LossFragmentError('Lack necessary fragment.')
        else:
            self.amide = Fragment(fragments['amide'])
            self.fragments.append(self.amide)

        if 'aryl' not in fragments:
            raise LossFragmentError('Lack necessary fragment.')
        else:
            self.aryl = Fragment(fragments['aryl'])
            self.fragments.append(self.aryl)

        if 'alcohol1' not in fragments and 'alcohol2' in fragments:
            raise LossFragmentError('Lack necessary fragment.')
        else:
            self.alcohol1 = Fragment(fragments['alcohol1'])
            self.fragments.append(self.alcohol1)

        if 'alcohol2' in fragments:
            self.alcohol2 = Fragment(fragments['alcohol2'])
            self.fragments.append(self.alcohol2)

        if 'modifier' in fragments:
            self.modifier = Fragment(fragments['modifier'])
            self.fragments.append(self.modifier)

        # make new_molecule object
        fragments_name_list = [fragment.name for fragment in self.fragments]
        new_mol_name = "_".join(fragments_name_list)
        self.new_molecule = Fragment(new_mol_name, new_molecule=True)

    def synthesize_attypes(self):
        new_attypes = []
        for fragment in self.fragments:
            new_attypes.extend(fragment.attypes[1:])
        new_attypes.insert(0, "")
        self.new_molecule.attypes = new_attypes

    # def synthesize_attypes(self):
    #     self.new_molecule.attypes = self.fragment1.attypes + self.fragment2.attypes[1:]

    def synthesize_bond(self):
        new_bonds = []
        for fragment in self.fragments:
            if not new_bonds:
                new_bonds = fragment.bond
                continue

            atom_max_idx = max(new_bonds)

            for bond in fragment.bond:
                new_bond = [i + atom_max_idx for i in bond]
                new_bonds.append(new_bond)

            adding_bond = self.get_adding_bond(fragment, atom_max_idx)
            new_bonds.append(adding_bond)

        self.new_molecule.bond = new_bonds

    def get_adding_bond(self, fragment: Fragment, atom_max_idx: int) -> List[List[int]]:

        bonds = defaultdict(list)
        adding_bond = []
        for i, bind_fragment in enumerate(fragment.bind_fragment):
            if bind_fragment:
                bond_type = sorted([self.get_fragment_type(fragment), bind_fragment])
                bond_type = tuple(bond_type)
                new_atom_idx = i + atom_max_idx
                bonds[bond_type].append(new_atom_idx)

        for bond_to_add in bonds.values():
            if len(bond_to_add) == 2:
                adding_bond.append(bond_to_add)

        return adding_bond


    # def synthesize_bond(self):
    #     fragment1_element_num = len(self.fragment1.attypes) - 1
    #     updated_list = []
    #     for bond in self.fragment2.bond:
    #         updated_list.append([i + fragment1_element_num for i in bond])
    #     self.new_molecule.bond = self.fragment1.bond + updated_list

    # def synthesize_free_atom(self):
    #     fragment1_element_num = len(self.fragment1.attypes) - 1
    #     adding_fragment_free_atom = [i + fragment1_element_num for i in self.fragment2.free_atom]
    #     self.new_molecule.free_atom = self.fragment1.free_atom + adding_fragment_free_atom

    # def synthesize_long_bond(self):
    #     fragment1_bond_num = len(self.fragment1.bond)
    #     adding_fragment_long_bond_idx = [i + fragment1_bond_num for i in self.fragment2.long_bond_idx]
    #     self.new_molecule.long_bond_idx = self.fragment1.long_bond_idx + adding_fragment_long_bond_idx

    @staticmethod
    def get_fragment_x_max(coord: List[List[float]]) -> float:
        """
        Return max x of fragment
        """
        coord = np.array(coord[1:])
        x_max = float(np.max(coord, axis=0)[0])
        return x_max

    @staticmethod
    def get_fragment_x_min(coord: List[List[float]]) -> float:
        """
        Return min x of fragment
        """
        coord = np.array(coord[1:])
        x_min = float(np.min(coord, axis=0)[0])
        return x_min

    def get_x_diff(self, coord1: List[List[float]], coord2: List[List[float]]) -> float:
        f1_max_x = self.get_fragment_x_max(coord1)
        f2_min_x = self.get_fragment_x_min(coord2)
        diff = f1_max_x - f2_min_x
        return diff

    def synthesize_coord(self):
        new_coord = []

        for fragment in self.fragments:
            if not new_coord:
                new_coord = fragment.coord
                continue
            diff = self.get_x_diff(new_coord, fragment.coord)
            for coord in fragment.coord:
                new_position = np.array(coord[1:]) + np.array([diff+1, 0, 0])
                new_coord.append(new_position)

            self.new_molecule.coord = new_coord

    # def synthesize_coord(self):
    #     diff = self.get_x_diff()
    #     coord_f1 = self.fragment1.coord
    #
    #     coord_f2 = np.array(self.fragment2.coord[1:])
    #     coord_f2 += np.array([diff + 1, 0, 0])
    #     coord_f2 = list(coord_f2)
    #     self.new_molecule.coord = coord_f1 + coord_f2

    def synthesize_bind_fragment(self):
        new_bind_fragment = []
        for fragment in self.fragments:
            new_bind_fragment.extend(fragment.bind_fragment[1:])
        new_bind_fragment.insert(0, "")
        self.new_molecule.attypes = new_bind_fragment

    # def add_bond(self):
    #     benzothiazole_amide = []
    #     amide_aryl = []
    #     aryl_alcohol1 = []
    #     aryl_alcohol2 = []
    #     benzothiazole_modifier = []
    #
    #     fragment_atom_num = 0
    #     fragment_atom_num_dict = {}
    #     for fragment in self.fragments:
    #         fragment_atom_num_dict.update(fragment=fragment_atom_num)
    #         fragment_atom_num += len(fragment)
    #
    #     for fragment in self.fragments:

    def get_fragment_type(self, fragment: Fragment) -> str:
        if fragment == self.benzothiazole:
            return 'benzothiazole'
        if fragment == self.amide:
            return 'amide'
        if fragment == self.aryl:
            return 'aryl'
        if fragment == self.alcohol1:
            return 'alcohol1'
        if fragment == self.alcohol2:
            return 'alcohol2'
        if fragment == self.modifier:
            return 'modifier'
