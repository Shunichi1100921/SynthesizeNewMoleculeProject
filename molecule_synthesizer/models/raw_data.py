import datetime
import os.path
import pathlib


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

        contents = [atom[0].upper() for atom in contents]
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
