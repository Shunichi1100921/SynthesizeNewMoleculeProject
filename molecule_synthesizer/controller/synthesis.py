"""Controller for synthesize fragment data and create necessary files."""
from molecule_synthesizer.models import file_data
from molecule_synthesizer.models.fragment import Fragment, Synthesis
from molecule_synthesizer.views import view


def synthesize_molecule_for_machine_learning() -> None:
    """Synthesize molecule data and create necessary files in 'Data/new_molecule_{DATE}{MONTH}_{YEAR}'
    The molecules to be synthesized are set in settings.py.

    ex) fragments = {'benzothiazole': 'F4', 'amide': 'F1', 'aryl': 'F5',
     'alcohol1': 'F2', 'alcohol2': 'F25', 'modifier': 'F13'}
    """
    fragment_sets = file_data.FragmentSet.load_data()

    for fragment_set in fragment_sets:
        synthesizer = Synthesis(fragment_set)
        synthesizer.remove_hydrogen()
        synthesizer.calculate_fragid()
        synthesizer.synthesize_attypes()
        synthesizer.synthesize_bond()
        synthesizer.synthesize_coord()
        # Save data
        new_mol = synthesizer.new_molecule
        create_file(new_mol)


def create_file(fragment: Fragment) -> None:
    """Create necessary files for machine learning indicated below.
        'attypes.dat'
        'bond.dat'
        'coord.dat'
        '{molecule_name}.pbd'
        '{molecule_name}.xyz'
        'fragid.dat'

    Args:
        fragment (Fragment): Fragment objects of new molecule.

    Returns:
        None
    """
    # contentの作成
    contents_creator = view.FileContentsCreator(fragment)
    attype_contents = contents_creator.get_attypes_contents()
    bond_contents = contents_creator.get_bond_contents()
    coord_contents = contents_creator.get_coord_contents()
    xyz_contents = contents_creator.get_xyz_file_contents()
    pdb_contents = contents_creator.get_pdb_file_contents()
    fragid_contents = contents_creator.get_fragid_contents()

    # Data Modelのインスタンス化の作成
    attype_model = file_data.AttypeData(fragment.name, new_molecule=True)
    bond_model = file_data.BondData(fragment.name, new_molecule=True)
    coord_model = file_data.CoordData(fragment.name, new_molecule=True)
    xyz_model = file_data.XYZFile(fragment.name, new_molecule=True)
    pdb_model = file_data.PDBFile(fragment.name, new_molecule=True)
    fragid_model = file_data.Fragid(fragment.name, new_molecule=True)

    # DataのSave
    attype_model.save_data(attype_contents)
    bond_model.save_data(bond_contents)
    coord_model.save_data(coord_contents)
    xyz_model.save_data(xyz_contents)
    pdb_model.save_data(pdb_contents)
    fragid_model.save_data(fragid_contents)
