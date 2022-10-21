from molecule_synthesizer.controller import chemical_reaction
from molecule_synthesizer.models.fragment import Fragment
from molecule_synthesizer.views import view


class LossFragmentError(Exception):
    pass


def synthesize_multiple_molecules() -> None:
    fragments_list = [
        ['F1', 'F2'],
        ['F1', 'F2', 'F10'],
        ['F3', 'F5', 'F14'],
        ['F6', 'F13', 'F25'],
        ['F12', 'F19', 'F32'],
        ['F15', 'F23', 'F36'],
    ]

    for fragments in fragments_list:
        synthesize_fragment(fragments)

    return None


def synthesize_fragment(fragments):
    # fragments = input('Input fragment type with space, like "F1 F2 F12"').split()
    new_mol = chemical_reaction.load_fragment_data(fragments.pop(0))

    while fragments:
        new_mol_removed_hydrogen = chemical_reaction.remove_hydrogen(new_mol)
        adding_mol = chemical_reaction.load_fragment_data(fragments.pop(0))
        adding_mol_removed_hydrogen = chemical_reaction.remove_hydrogen(adding_mol)
        new_mol = chemical_reaction.synthesize_two_fragment(new_mol_removed_hydrogen, adding_mol_removed_hydrogen)

    create_file(new_mol)
    return new_mol


def synthesize_new_mol_for_machine_learning(modifier: str, benzothiazole: str, amide: str,
                                                aryl: str, alcohol1: str, alcohol2: str) -> Fragment:

    new_mol = Fragment(f'{modifier}_{benzothiazole}_{amide}_{aryl}_{alcohol1}_{alcohol2}')

    modifier = Fragment(modifier)
    benzothiazole = Fragment(benzothiazole)
    amide = Fragment(amide)
    alcohol1 = Fragment(alcohol1)
    alcohol2 = Fragment(alcohol2)
    aryl = Fragment(aryl)

    if not benzothiazole or not amide or not aryl:
        raise LossFragmentError('New molecule all needs Benzothiazole, Amide, and Aryl group.')

    benzothiazole.remove_hydrogen('amide')
    amide.remove_hydrogen('benzothiazole')
    new_mol.synthesize(benzothiazole)
    new_mol.synthesize(amide)

    new_mol.remove_hydrogen('aryl')
    aryl.remove_hydrogen('amide')
    new_mol.synthesize(aryl)

    if alcohol1:
        new_mol.remove_hydrogen('alcohol1')
        alcohol1.remove_hydrogen('aryl')
        new_mol.synthesize(alcohol1)

    if alcohol2:
        new_mol.remove_hydrogen('alcohol2')
        alcohol2.remove_hydrogen('aryl')
        new_mol.synthesize(alcohol2)

    if modifier:
        new_mol.remove_hydrogen('modifier')
        modifier.remove_hydrogen('benzothiazole')
        new_mol.synthesize(modifier)
    return new_mol


def create_file(fragment_data):
    # contentの作成
    contents_creator = view.FileContentsCreator(fragment_data)
    attype_contents = contents_creator.get_attypes_contents()
    bond_contents = contents_creator.get_bond_contents()
    coord_contents = contents_creator.get_coord_contents()
    xyz_contents = contents_creator.get_xyz_file_contents()
    pdb_contents = contents_creator.get_pdb_file_contents()

    # Data Modelのインスタンス化の作成
    attype_model = fragment_data.AttypeData(fragment_data['name'], new_molecule=True)
    bond_model = fragment_data.BondData(fragment_data['name'], new_molecule=True)
    coord_model = fragment_data.CoordData(fragment_data['name'], new_molecule=True)
    xyz_model = fragment_data.XYZFile(fragment_data['name'], new_molecule=True)
    pdb_model = fragment_data.PDBFile(fragment_data['name'], new_molecule=True)

    # DataのSave
    attype_model.save_data(attype_contents)
    bond_model.save_data(bond_contents)
    coord_model.save_data(coord_contents)
    xyz_model.save_data(xyz_contents)
    pdb_model.save_data(pdb_contents)
