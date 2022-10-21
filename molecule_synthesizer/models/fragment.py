from typing import List

import numpy as np

from molecule_synthesizer.models import fragment_data, chemical_synthesis


class FragmentModel(object):
    def __init__(self, fragment_name: str = None) -> None:
        # Make instance of data model.
        self.name = fragment_name

        self.attypes = fragment_data.AttypeData(self.name)
        self.bond = fragment_data.BondData(self.name)
        self.coord = fragment_data.CoordData(self.name)
        self.bind_fragment = fragment_data.BindFragmentData(self.name)
        self.xyz = fragment_data.XYZFile(self.name)
        self.pdb = fragment_data.PDBFile(self.name)


    def remove_hydrogen(self, bind_fragment: str) -> None:
        """Remove Hydrogen which is seemed to be added.
        Process:
            1st: Remove bond from bond data.
            2nd: Remove a hydrogen data from attypes data.
            3rd: Arrange bond index of bond data.
            4th: Remove a coord data of removed hydrogen from coord data.
            5th: Arrange a bond index of free-atom data.

        :param fragment_data: dict: fragment name such as 'F1', 'F2', ..., 'F57'.
        :return: new_fragment_data(dict): Edited fragment data that shows attypes, bond, long bond, free atom, and coord.
        """
        remover = chemical_synthesis.RemoveHydrogen(self)
        remover.cut_bond()
        remover.remove_hydrogen_from_attype()
        remover.arrange_bond_idx()
        remover.remove_hydrogen_from_coord()
        remover.arrange_free_atom_idx()


class Fragment(FragmentModel):
     def __init__(self, fragment_name: str = None) -> None:



        super().__init__(fragment_name)
        self.attypes_data = self.attypes.load_data()
        self.bond_data = self.bond.load_data()
        self.coord_data = self.coord.load_data()
        self.bind_fragment_data = self.bind_fragment.load_data()
        self.free_atom = None


class NewMolecule(FragmentModel):

    def __init__(self, fragment_name: str = None) -> None:
        self.name = None
        self.attypes_data = None
        self.bond_data = None
        self.coord_data = None
        self.bind_fragment_data = None
        self.free_atom = None
        super().__init__(fragment_name)

    def synthesize(self, fragment: Fragment) -> None:
        if not self.name:
            self.name = fragment.name
            self.attypes_data = fragment.attypes_data
            self.bond_data = fragment.bond_data
            self.coord_data = fragment.coord_data
            self.bind_fragment_data = fragment.bind_fragment_data
            self.free_atom = fragment.free_atom
        else:
            self.synthesize_name(fragment)
            self.synthesize_attypes(fragment)
            self.arrange_atom_index_of_bond(fragment)
            self.arrange_atom_index_of_free_atom(fragment)
            self.synthesize_bond(fragment)
            self.synthesize_free_atom(fragment)
            self.arrange_coord_to_place_fragment(fragment)

    def synthesize_name(self, fragment: Fragment) -> None:
        self.name = f'{self.name}_{fragment.name}'

    def synthesize_attypes(self, fragment: Fragment) -> None:
        self.attypes_data += fragment.attypes_data[1:]

    def arrange_atom_index_of_bond(self, fragment: Fragment) -> None:
        element_num = len(self.attypes_data) - 1
        updated_list = []
        for bond in fragment.bond_data:
            updated_list.append((i + element_num for i in bond))

        fragment.bond_data = updated_list

    def arrange_atom_index_of_free_atom(self, fragment: Fragment) -> None:
        element_num = len(self.attypes_data) - 1
        fragment.free_atom = [i + element_num for i in fragment.free_atom]

    def synthesize_bond(self, fragment: Fragment):
        self.bond_data += fragment.bond_data

    def synthesize_free_atom(self, fragment: Fragment):
        self.free_atom += fragment.free_atom

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

    def get_x_diff(self, fragment: Fragment) -> float:
        new_mol_max_x = self.get_fragment_x_min_max(self.coord_data)[1]
        fragment_min_x = self.get_fragment_x_min_max(fragment.coord_data)[0]
        if new_mol_max_x < fragment_min_x:
            return -1
        diff = new_mol_max_x - fragment_min_x
        return diff

    def arrange_coord_to_place_fragment(self, fragment: Fragment):
        diff = self.get_x_diff(fragment)
        fragment_coord = np.array(fragment.coord_data[1:])
        fragment_coord += np.array([diff+1, 0, 0])
        fragment_coord = fragment_coord.tolist()
        self.coord_data += fragment_coord

    def add_bond(self):
        free_atom_1 = self.free_atom.pop(0)
        free_atom_2 = self.free_atom.pop(0)
        new_bond = [free_atom_1, free_atom_2]
        self.bond_data.append(new_bond)











