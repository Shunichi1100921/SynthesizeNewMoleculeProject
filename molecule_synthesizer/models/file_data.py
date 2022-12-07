from __future__ import annotations

import datetime
import os.path
import pathlib

import settings


class FragmentDataModel(object):
    def __init__(self, fragment_type, new_molecule=False):
        self.fragment_type = fragment_type
        self.data_dir_path = self.get_dir_path(fragment_type, new_molecule)
        self.file_path = self.get_file_path()

    @staticmethod
    def get_dir_path(fragment_type, new_molecule):
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

    def save_data(self, contents):
        with open(self.file_path, 'w') as f:
            f.writelines(contents)


class AttypeData(FragmentDataModel):
    def get_file_path(self):
        file_path = os.path.join(self.data_dir_path, 'attypes.dat')
        return file_path

    def load_data(self):

        with open(self.file_path, 'r') as f:
            contents = f.readlines()

        contents = [atom.split('_')[0] for atom in contents]
        contents.insert(0, '')
        return contents


class BondData(FragmentDataModel):
    def get_file_path(self):
        file_path = os.path.join(self.data_dir_path, 'bond.dat')
        return file_path

    def load_data(self):
        with open(self.file_path, 'r') as f:
            contents = f.readlines()
        data = [bond.split() for bond in contents]

        # long bondをbondのindexとして保存する
        long_bond_list = []
        for i, bond in enumerate(data):
            if len(bond) >= 3:
                long_bond_list.append(i)

        # 数字のリストに変換する ex) [['1', '2', '*****'], ['2', '12']] -> [[1, 2], [2, 12]]
        data = [bond[:2] for bond in data]
        data = [list(map(int, bond)) for bond in data]
        return long_bond_list, data


class CoordData(FragmentDataModel):
    def get_file_path(self):
        file_path = os.path.join(self.data_dir_path, 'coord.dat')
        return file_path

    def load_data(self):
        with open(self.file_path, 'r') as f:
            contents = f.readlines()

        contents = [list(map(float, coord.split())) for coord in contents]
        contents.insert(0, [])
        return contents


class BindFragment(FragmentDataModel):
    def get_file_path(self):
        file_path = os.path.join(self.data_dir_path, 'bind_fragment.dat')
        return file_path

    def load_data(self):
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
    def get_file_path(self):
        frag_num = self.fragment_type[1:]
        file_path = os.path.join(self.data_dir_path, f'F{frag_num}.xyz')
        return file_path


class PDBFile(FragmentDataModel):
    def get_file_path(self):
        frag_type = self.fragment_type
        file_path = os.path.join(self.data_dir_path, f'{frag_type}.pdb')
        return file_path


class Fragid(FragmentDataModel):
    def get_file_path(self):
        file_path = os.path.join(self.data_dir_path, 'fragid.dat')
        return file_path


class FragmentSet(object):
    @classmethod
    def load_data(cls) -> list[dict[str, str | None]]:
        """Load Fragment sets
        Returns:
            fragment_sets(list): A list of dictionary that create one new molecule.
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
            for k, v in fragments.items():
                if not type(v) == list:
                    fragments[k] = [v]

        for benzothiazole in fragments['benzothiazole']:
            for amide in fragments['amide']:
                for aryl in fragments['aryl']:
                    for alcohol1 in fragments['alcohol1']:
                        for alcohol2 in fragments['alcohol2']:
                            for modifier in fragments['modifier']:
                                if aryl == 'F55' or aryl == 'F56':
                                    alcohol1, alcohol2 = None, None
                                if benzothiazole == 'F15' or benzothiazole == 'F57':
                                    modifier = None
                                if alcohol1 == None:
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
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        file_path = os.path.join(base_dir, 'Data', 'fragment_classification.dat')
        fragments = {}
        with open(file_path, 'r') as f:
            contents = f.read()

        contents = contents.split('#')
        contents = [l.rstrip(os.linesep).split('\n') for l in contents]
        for l in contents:
            fragments.update({l[0].lower().lstrip(): l})
        fragments.update({'alcohol1': fragments['alcohol'], 'alcohol2': fragments['alcohol']})

        fragments['alcohol1'].append(None)
        fragments['alcohol2'].append(None)
        fragments['modifier'].append(None)

        return fragments
