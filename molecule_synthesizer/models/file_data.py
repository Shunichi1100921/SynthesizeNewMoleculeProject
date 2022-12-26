"""Open and process and save data."""
from __future__ import annotations

import datetime
import os.path
import pathlib

import settings


class FragmentDataModel(object):
    """Processing data of file.
    """
    def __init__(self, fragment_type, new_molecule=False) -> None:
        self.fragment_type = fragment_type
        self.data_dir_path = self.get_dir_path(fragment_type, new_molecule)
        self.file_path = self.get_file_path()

    @staticmethod
    def get_dir_path(fragment_type, new_molecule) -> str:
        """Get directory path of data of fragment.

        Args:
            fragment_type (str) : 'F1', 'F2', 'F3', ...
            new_molecule (bool) : If molecule does already exist, False, otherwise, True.
        Return:
            data_dir_path: '.../F1', '.../F2', ...
        """

        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        if new_molecule:
            today = datetime.date.today().strftime('%d%b_%Y')
            data_base_dir_path = os.path.join(base_dir, 'Data', f'new_molecule_{today}')
            if not os.path.exists(data_base_dir_path):
                pathlib.Path(data_base_dir_path).mkdir()
        else:
            data_base_dir_path = os.path.join(base_dir, 'Data', 'fragment_data_for_shape_25Aug_2022')

        data_dir_path = os.path.join(data_base_dir_path, fragment_type)

        if new_molecule and not os.path.exists(data_dir_path):
            pathlib.Path(data_dir_path).mkdir()

        return data_dir_path

    def get_file_path(self):
        pass

    def load_data(self):
        pass

    def save_data(self, contents) -> None:
        """Save data.  The type of data should be list."""
        with open(self.file_path, 'w') as f:
            f.writelines(contents)


class AttypeData(FragmentDataModel):
    """Processing Attype data which means 'atom type'.  In each molecule, the type of each atom is written."""
    def get_file_path(self) -> str:
        file_path = os.path.join(self.data_dir_path, 'attypes.dat')
        return file_path

    def load_data(self) -> list[str]:
        """Open file and process the index of the list to match the index of the atom.

        Returns:
            contents: List of atom type, i.g., ['C', 'C', 'N', 'H'].
            """

        with open(self.file_path, 'r') as f:
            contents = f.readlines()

        contents = [atom.split('_')[0] for atom in contents]
        contents.insert(0, '')
        return contents


class BondData(FragmentDataModel):
    """Processing bond data.  Each atom of each fragment is assigned an index for each fragment,
     and bonds are represented by that index."""
    def get_file_path(self) -> str:
        file_path = os.path.join(self.data_dir_path, 'bond.dat')
        return file_path

    def load_data(self) -> tuple[list[int], list[list[int]]]:
        """Open file and store index of atoms of each bond in list.

        Returns:
            long_bond_list (list): Long bonds, i.e., the bonds added later,
                                    is stored in the list as the index of bond list.
            data (list): Return the list containing the list which store the index of the atoms of each bond.

        """
        with open(self.file_path, 'r') as f:
            contents = f.readlines()
        data = [bond.split() for bond in contents]

        # Save bond index which seems to be added later.
        long_bond_list = []
        for i, bond in enumerate(data):
            if len(bond) >= 3:
                long_bond_list.append(i)

        # Change type of data from strings to integer and remove
        # unneeded parts. ex) [['1', '2', '*****'], ['2', '12']] -> [[1, 2], [2, 12]]
        data = [bond[:2] for bond in data]
        data = [list(map(int, bond)) for bond in data]
        return long_bond_list, data


class CoordData(FragmentDataModel):
    """Coord data represents the position of each atom."""
    def get_file_path(self) -> str:
        file_path = os.path.join(self.data_dir_path, 'coord.dat')
        return file_path

    def load_data(self) -> list[list[float]]:
        """Open file and process the index of the list to match the index of the atom.

        Returns:
            contents (list): Float data in which the position of each atom is expressed in xyz coordinates.

        """
        with open(self.file_path, 'r') as f:
            contents = f.readlines()

        contents = [list(map(float, coord.split())) for coord in contents]
        contents.insert(0, [])
        return contents


class BindFragment(FragmentDataModel):
    """Bind Fragment shows which type of fragment that each atom binds to."""
    def get_file_path(self) -> str:
        file_path = os.path.join(self.data_dir_path, 'bind_fragment.dat')
        return file_path

    def load_data(self) -> list[str]:
        """Open file, lowercase each data, and align the indexes of the list and atoms.

        Returns:
            contents(list): Fragment types are listed.  It contains "", which means there won't bind other fragments.

        """
        try:
            with open(self.file_path, 'r') as f:
                contents = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f'There is not bind_fragment.dat file in {self.fragment_type}\n'
                                    f'The file was created manually, so there may be omissions.'
                                    f'Please create bind_fragment.dat file in {self.fragment_type}')

        contents = [item.rstrip(os.linesep) for item in contents]
        contents = [item.lower() for item in contents]
        contents.insert(0, '')
        return contents


class XYZFile(FragmentDataModel):
    """XYZ File is one of the output format for expressing molecular."""
    def get_file_path(self) -> str:
        frag_num = self.fragment_type[1:]
        file_path = os.path.join(self.data_dir_path, f'F{frag_num}.xyz')
        return file_path


class PDBFile(FragmentDataModel):
    """PDB File (means, Protein DataBase) is another output format for expressing molecular.
      This is more popular than XYZ."""

    def get_file_path(self) -> str:
        frag_type = self.fragment_type
        file_path = os.path.join(self.data_dir_path, f'{frag_type}.pdb')
        return file_path


class Fragid(FragmentDataModel):
    """Fragid is necessary for machine learning.  It shows which fragment each atom derived from."""
    def get_file_path(self):
        file_path = os.path.join(self.data_dir_path, 'fragid.dat')
        return file_path


class FragmentSet(object):
    @classmethod
    def load_data(cls) -> list[dict[str, str | None]]:
        """Load Fragment sets from settings.py
        Returns:
            fragment_sets(list): A list of dictionary. Each dictionary have fragment name to create one new molecule.
        """
        fragment_sets = []
        if settings.All_Fragment:
            fragments = cls.load_all_fragments()
        else:
            fragments = {
                'benzothiazole': settings.benzothiazole,
                'amide': settings.amide,
                'aryl': settings.aryl,
                'alcohol1': settings.alcohol1,
                'alcohol2': settings.alcohol2,
                'modifier': settings.modifier
                         }
            # Align the type of value in the dictionary with list.
            for k, v in fragments.items():
                if not type(v) == list:
                    fragments[k] = [v]

        # Search all patterns in a linear search.
        for benzothiazole in fragments['benzothiazole']:
            for amide in fragments['amide']:
                for aryl in fragments['aryl']:
                    for alcohol1 in fragments['alcohol1']:
                        for alcohol2 in fragments['alcohol2']:
                            for modifier in fragments['modifier']:
                                if aryl in ['F55', 'F56']:
                                    alcohol1, alcohol2 = None, None
                                if benzothiazole in ['F15', 'F57']:
                                    modifier = None
                                if alcohol1 is None:
                                    alcohol2 = None
                                fragment_set = {
                                    'benzothiazole': benzothiazole,
                                    'amide': amide,
                                    'aryl': aryl,
                                    'alcohol1': alcohol1,
                                    'alcohol2': alcohol2,
                                    'modifier': modifier,
                                }
                                fragment_sets.append(fragment_set)

        return fragment_sets

    @classmethod
    def load_all_fragments(cls) -> dict[str, list[str | None]]:
        """Load all fragments from 'fragment_classification.dat'

        Returns:
            fragments (dict): For each fragment type, the value is a list of fragment names.

        """
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        file_path = os.path.join(base_dir, 'Data', 'fragment_classification.dat')
        fragments = {}
        with open(file_path, 'r') as f:
            contents = f.read()

        contents = contents.split('#')
        contents = [fragment_name.rstrip(os.linesep).split('\n') for fragment_name in contents]
        for line in contents:
            fragments.update({line[0].lower().lstrip(): line[1:]})
        fragments.update({'alcohol1': fragments['alcohol'], 'alcohol2': fragments['alcohol']})

        fragments['alcohol1'].append(None)
        fragments['alcohol2'].append(None)
        fragments['modifier'].append(None)

        return fragments
