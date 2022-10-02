from molecule_synthesizer.controller import chemical_reaction
from molecule_synthesizer.models import raw_data
from molecule_synthesizer.views import view


def synthesize_fragment():
    # fragments = input('Input fragment type with space, like "F1 F2 F12"').split()
    fragments = ['F12', 'F19', 'F32']
    new_mol = chemical_reaction.load_fragment_data(fragments.pop(0))

    while fragments:
        new_mol_removed_hydrogen = chemical_reaction.remove_hydrogen(new_mol)
        adding_mol = chemical_reaction.load_fragment_data(fragments.pop(0))
        adding_mol_removed_hydrogen = chemical_reaction.remove_hydrogen(adding_mol)
        new_mol = chemical_reaction.synthesize_two_fragment(new_mol_removed_hydrogen, adding_mol_removed_hydrogen)

    create_file(new_mol)
    return new_mol

def create_file(fragment_data):
    contents_creator = view.FileContentsCreator(fragment_data)
    attype_model = raw_data.AttypeData(fragment_data['name'], new_molecule=True)
    bond_model = raw_data.BondData(fragment_data['name'], new_molecule=True)
    coord_model = raw_data.CoordData(fragment_data['name'], new_molecule=True)
    xyz_model = raw_data.XYZFile(fragment_data['name'], new_molecule=True)
    pdb_model = raw_data.PDBFile(fragment_data['name'], new_molecule=True)

    attype_contents = contents_creator.get_attypes_contents()
    bond_contents = contents_creator.get_bond_contents()
    coord_contents = contents_creator.get_coord_contents()
    xyz_contents = contents_creator.get_xyz_file_contents()
    pdb_contents = contents_creator.get_pdb_file_contents()

    attype_model.save_data(attype_contents)
    bond_model.save_data(bond_contents)
    coord_model.save_data(coord_contents)
    xyz_model.save_data(xyz_contents)
    pdb_model.save_data(pdb_contents)
