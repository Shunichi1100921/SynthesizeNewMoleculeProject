"""Chemical reaction edit 'attypes', 'bond', 'coord', and 'free_atom' data."""
from __future__ import annotations

from collections import defaultdict
from typing import List, Dict, Tuple
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


Bond = List[int]
Coord = List[float]


class Fragment(object):
    """Each fragment is treated as a unit of this class."""

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
        self.fragid = []

    def load_data(self) -> None:
        """Load preprocessed data."""
        self.attypes = self.attypes_obj.load_data()
        self.long_bond_idx, self.bond = self.bond_obj.load_data()
        self.coord = self.coord_obj.load_data()
        self.bind_fragment = self.bind_fragment_obj.load_data()

    def remove_hydrogen(self, bind_fragment: str) -> None:
        """Remove Hydrogen which for binding another bond.
        Process:
            1st: Remove a bond with hydrogen which is long and seemed to be added from bond data.
            2nd: Remove a hydrogen data from attypes data.
            3rd: Arrange bond index of bond data.
            4th: Remove a coord data of removed hydrogen from coord data.
            5th: Arrange a bond index of free-atom data.

        Args:
            bind_fragment (str): Fragment type the fragment will bind with.

        """

        remover = RemoveHydrogen(self, bind_fragment)
        remover.remove_bond()


class RemoveHydrogen(object):
    """Removing hydrogen process Fragment object data."""

    def __init__(self, fragment: Fragment, bind_fragment: str) -> None:
        self.fragment = fragment
        self.bind_fragment = bind_fragment
        self.bond_idx_to_remove = self.find_bond_idx_to_remove()
        self.atom_idx_by_attype = self.find_atom_idx_by_attype()

    def find_bond_idx_to_remove(self) -> int:
        """Return one bond to cut.

        Return:
            bond_idx: The index of bond to remove in bonds list .
        Raises:
            BindFragmentError: Raise this error when the fragment don't have appropriate bind_fragment data.
        """

        for bond_idx in self.fragment.long_bond_idx:
            bond = self.fragment.bond[bond_idx]
            for atom_idx in bond:
                if atom_idx >= len(self.fragment.bind_fragment):
                    continue
                if self.fragment.bind_fragment[atom_idx] == self.bind_fragment:
                    return bond_idx

        raise BindFragmentError(f"There is not appropriate fragment.  Fragments are {self.fragment.name}")

    def find_atom_idx_by_attype(self) -> dict[str, int]:
        """Returns which of the bonds to remove is hydrogen and will bind with 'bind_fragment'.
        Return:
            attype_in_bond (dict): Returns a dictionary with key as hydrogen and {bind_fragment},
            i.e., {'hydrogen': atom_idx, {bind_fragment}: atom_idx}
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
        """Remove bond and hydrogen and arrange index.

        """
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
        6th: add bond <- synthesize with "synthesize bond"
            1st: synthesize atom idx of bond
            2nd: add bond
            3rd: remove free atom
        7th: synthesize bind_fragment <- Don't need?
    """

    def __init__(self, fragments: dict[str, str]) -> None:
        """Create Fragment dataset.
        Args:
            fragments: Dictionary of fragment part and fragment name.
                ex)  {'benzothiazole': 'F4', 'amide': 'F1', 'aryl': 'F5',
                'alcohol1': 'F2', 'alcohol2': 'F25', 'modifier': 'F13'}
        """
        self.fragments = []

        # Make fragment object.

        if fragments['modifier']:
            self.modifier = Fragment(fragments['modifier'])
            self.fragments.append(self.modifier)
        else:
            self.modifier = None

        if not fragments['benzothiazole']:
            raise LossFragmentError('Lack necessary fragment.')
        else:
            self.benzothiazole = Fragment(fragments['benzothiazole'])
            self.fragments.append(self.benzothiazole)

        if not fragments['amide']:
            raise LossFragmentError('Lack necessary fragment.')
        else:
            self.amide = Fragment(fragments['amide'])
            self.fragments.append(self.amide)

        if not fragments['aryl']:
            raise LossFragmentError('Lack necessary fragment.')
        else:
            self.aryl = Fragment(fragments['aryl'])
            self.fragments.append(self.aryl)

        if not fragments['alcohol1']:
            if fragments['alcohol2']:
                raise LossFragmentError('Lack necessary fragment.')
            else:
                self.alcohol1 = None
        else:
            self.alcohol1 = Fragment(fragments['alcohol1'])
            self.fragments.append(self.alcohol1)

        if fragments['alcohol2']:
            self.alcohol2 = Fragment(fragments['alcohol2'])
            self.fragments.append(self.alcohol2)
        else:
            self.alcohol2 = None

        # make new_molecule object
        fragments_name_list = [fragment.name for fragment in self.fragments]
        new_mol_name = "_".join(fragments_name_list)
        self.new_molecule = Fragment(new_mol_name, new_molecule=True)

    def remove_hydrogen(self) -> None:
        """Call remove_hydrogen method of fragment object with appropriate args."""
        self.benzothiazole.remove_hydrogen('amide')
        self.amide.remove_hydrogen('benzothiazole')
        self.amide.remove_hydrogen('aryl')
        self.aryl.remove_hydrogen('amide')

        if self.alcohol1:
            self.aryl.remove_hydrogen('alcohol1')
            self.alcohol1.remove_hydrogen('aryl')
        if self.alcohol2:
            self.aryl.remove_hydrogen('alcohol2')
            self.alcohol2.remove_hydrogen('aryl')
        if self.modifier:
            self.modifier.remove_hydrogen('benzothiazole')
            self.benzothiazole.remove_hydrogen('modifier')

    def calculate_fragid(self) -> None:
        """Count the number of atoms in each fragment.
        Notes: It should be performed AFTER THE REMOVE_HYDROGEN METHOD EXECUTED,
         since it is the number of atoms after removing hydrogen that is important.
        """
        fragid = []
        for fragment in self.fragments:
            atom_num_of_fragment = len(fragment.attypes) - 1
            for i in range(atom_num_of_fragment):
                fragid.append(fragment.name)
        self.new_molecule.fragid = fragid

    def synthesize_attypes(self) -> None:
        """Synthesize attypes extending list.
        """

        new_attypes = []
        for fragment in self.fragments:
            new_attypes.extend(fragment.attypes[1:])
        new_attypes.insert(0, "")
        self.new_molecule.attypes = new_attypes

    def synthesize_bond(self):
        """Synthesize all bond data updating atom index and add new bond.
        First add the bond that was originally in each fragment, then add the new bond.
        """

        new_bonds = []
        created_bonds = defaultdict(list)
        for fragment in self.fragments:
            if new_bonds:
                bond = np.array(new_bonds)
                atom_max_idx = np.max(bond)
            else:
                atom_max_idx = 0

            # Add existing bonds.
            for bond in fragment.bond:
                updated_bond = [i + atom_max_idx for i in bond]
                new_bonds.append(updated_bond)

            # Update created bond.
            created_bonds = self.get_created_bond(fragment, atom_max_idx, created_bonds)

        # Add created bond.
        for bond_to_add in created_bonds.values():
            if len(bond_to_add) == 2:
                new_bonds.append(bond_to_add)

        self.new_molecule.bond = new_bonds

    def get_created_bond(self, fragment: Fragment, atom_max_idx: int, created_bonds: dict) -> dict[tuple[str], Bond]:
        """Edit dictionary of bonds created.

        Args:
            fragment (Fragment): About what fragment.
            atom_max_idx (int): Maximum index value of atoms to which bonds have already been added.
            Use this value for calculating new atom index of atom of new bond.
            created_bonds (dict): Key are Fragment types and value is the index of each atom.

        Returns:
            created_bond: Refer to args.

        """

        fragment_type = self.get_fragment_type(fragment)
        for i, bind_fragment in enumerate(fragment.bind_fragment):
            if bind_fragment:
                bond_fragment_type = sorted([fragment_type, bind_fragment])
                bond_fragment_type = tuple(bond_fragment_type)
                new_atom_idx = i + atom_max_idx
                created_bonds[bond_fragment_type].append(new_atom_idx)
        return created_bonds

    @staticmethod
    def get_fragment_x_max(coord: List[List[float]]) -> float:
        """Return max x value of the fragment.

        Args:
            coord: Positional data of atom.

        Returns:
            x_max: Max x-value of the fragment.

        """
        coord = np.array(coord[1:])
        x_max = float(np.max(coord, axis=0)[0])
        return x_max

    @staticmethod
    def get_fragment_x_min(coord: List[List[float]]) -> float:
        """Return min x value of the fragment.

        Args:
            coord: Positional data of atom.

        Returns:
            x_min: min x-value of the fragment.

        """
        coord = np.array(coord[1:])
        x_min = float(np.min(coord, axis=0)[0])
        return x_min

    def get_x_diff(self, coord1: List[List[float]], coord2: List[List[float]]) -> float:
        """Calculate the distance of the x-coordinate of atoms.

        Args:
            coord1: The positional data of one atom.
            coord2: The positional data of the other atom.

        Returns:
            diff (float): The distance of the x-cooridnate of atoms.

        """
        f1_max_x = self.get_fragment_x_max(coord1)
        f2_min_x = self.get_fragment_x_min(coord2)
        diff = f1_max_x - f2_min_x
        return diff

    def synthesize_coord(self) -> None:
        """Synthesize coord data in order not to overlap each atom.

        """
        new_coord = []

        for fragment in self.fragments:
            if not new_coord:
                new_coord = fragment.coord
                continue
            diff = self.get_x_diff(new_coord, fragment.coord)
            # for coord in fragment.coord:
            #     new_position = np.array(coord[1:]) + np.array([diff+1, 0, 0])
            new_position = np.array(fragment.coord[1:]) + np.array([diff+1, 0, 0])
            new_coord.extend(new_position)

            self.new_molecule.coord = new_coord

    def synthesize_bind_fragment(self) -> None:
        """Synthesize bind_fragment data arranging atom index.

        """
        new_bind_fragment = []
        for fragment in self.fragments:
            new_bind_fragment.extend(fragment.bind_fragment[1:])
        new_bind_fragment.insert(0, "")
        self.new_molecule.attypes = new_bind_fragment

    def get_fragment_type(self, fragment: Fragment) -> str:
        """Utility function for getting fragment type name as strings.

        Args:
            fragment(Fragment): Fragment object

        Returns:
            fragment_type(str): Fragment type in a new molecule.
        """
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
