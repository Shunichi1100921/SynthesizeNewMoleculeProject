from typing import Dict

from molecule_synthesizer.models import file_data
from molecule_synthesizer.models.fragment import Fragment, Synthesis
from molecule_synthesizer.views import view

# def synthesize_multiple_molecules() -> None:
#     fragments_list = [
#         ['F1', 'F2'],
#         ['F1', 'F2', 'F10'],
#         ['F3', 'F5', 'F14'],
#         ['F6', 'F13', 'F25'],
#         ['F12', 'F19', 'F32'],
#         ['F15', 'F23', 'F36'],
    # ]
    #
    # for fragments in fragments_list:
    #     synthesize_fragment(fragments)
    #
    # return None


# def synthesize_fragment(fragments):
#     fragments = input('Input fragment type with space, like "F1 F2 F12"').split()
    # new_mol = Fragment(fragments.pop(0))
    # new_mol.load_data()
    #
    # while fragments:
    #     new_mol.remove_hydrogen()
    #     adding_mol = Fragment(fragments.pop(0))
    #     adding_mol.load_data()
    #     adding_mol.remove_hydrogen()
    #     new_mol = chemical_reaction.synthesize_two_fragment(new_mol, adding_mol)
    #
    # create_file(new_mol)
    # return new_mol


def synthesize_molecule_for_machine_learning(fragments: Dict[str, str]) -> None:
    """
    ex) fragments = {'benzothiazole': 'F4', 'amide': 'F1', 'aryl': 'F5',
     'alcohol1': 'F2', 'alcohol2': 'F25', 'modifier': 'F13'}
    """

    # Fragmentの生成、removehydrogen、synthesizeまで全部synthesizerが行う
    # new_moleculeの情報は保存する
    synthesizer = Synthesis(fragments)
    synthesizer.synthesize('benzothiazole', 'amide')
    synthesizer.synthesize('amide', 'aryl')
    synthesizer.synthesize('aryl', 'alcohol1')
    synthesizer.synthesize('aryl', 'alcohol2')
    synthesizer.synthesize('benzothiazole', 'modifier')

    new_mol = synthesizer.new_molecule

    create_file(new_mol)

    return None

def create_file(fragment: Fragment) -> None:
    # contentの作成
    contents_creator = view.FileContentsCreator(fragment)
    attype_contents = contents_creator.get_attypes_contents()
    bond_contents = contents_creator.get_bond_contents()
    coord_contents = contents_creator.get_coord_contents()
    xyz_contents = contents_creator.get_xyz_file_contents()
    pdb_contents = contents_creator.get_pdb_file_contents()

    # Data Modelのインスタンス化の作成
    attype_model = file_data.AttypeData(fragment.name, new_molecule=True)
    bond_model = file_data.BondData(fragment.name, new_molecule=True)
    coord_model = file_data.CoordData(fragment.name, new_molecule=True)
    xyz_model = file_data.XYZFile(fragment.name, new_molecule=True)
    pdb_model = file_data.PDBFile(fragment.name, new_molecule=True)

    # DataのSave
    attype_model.save_data(attype_contents)
    bond_model.save_data(bond_contents)
    coord_model.save_data(coord_contents)
    xyz_model.save_data(xyz_contents)
    pdb_model.save_data(pdb_contents)
