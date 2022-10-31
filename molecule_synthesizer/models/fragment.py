"""Chemical reaction edit 'attypes', 'bond', 'coord', and 'free_atom' data."""
from typing import List, Dict
import numpy as np

from molecule_synthesizer.models import file_data


class BondError(Exception):
    pass

class LongBondError(Exception):
    pass

class BindFragmentError(Exception):
    pass


class Fragment(object):

    def __init__(self, name:str = None, new_molecule: bool = False) -> None:
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
            self.free_atom = []
            self.bind_fragment = []
        else:
            self.load_data()

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
