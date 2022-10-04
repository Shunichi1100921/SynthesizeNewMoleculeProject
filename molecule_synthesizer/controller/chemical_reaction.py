"""Data: atom types, bond, coord, long-bond, free-atom"""
from molecule_synthesizer.models import chemical_synthesis
from molecule_synthesizer.models import raw_data


def load_fragment_data(fragment):
    fragment_data = {}
    attypes = raw_data.AttypeData(fragment)
    bond = raw_data.BondData(fragment)
    coord = raw_data.CoordData(fragment)

    fragment_data['name'] = fragment
    fragment_data['attypes'] = attypes.load_data()
    fragment_data['long_bond'], fragment_data['bond'] = bond.load_data()
    fragment_data['coord'] = coord.load_data()
    fragment_data['free_atom'] = []
    return fragment_data


def remove_hydrogen(fragment_data):
    """Remove Hydrogen which is seemed to be added.
    Process:
        1st: Remove a bond with hydrogen which is long and seemed to be added from bond data.
        2nd: Remove a hydrogen data from attypes data.
        3rd: Arrange bond index of bond data.
        4th: Remove a coord data of removed hydrogen from coord data.
        5th: Arrange a bond index of free-atom data.

    :param fragment_data: dict: fragment name such as 'F1', 'F2', ..., 'F57'.
    :return: new_fragment_data(dict): Edited fragment data that shows attypes, bond, long bond, free atom, and coord.
    """
    remover = chemical_synthesis.RemoveHydrogen(fragment_data)
    remover.cut_bond()
    remover.remove_hydrogen_from_attype()
    remover.arrange_bond_idx()
    remover.remove_hydrogen_from_coord()
    remover.arrange_free_atom_idx()
    new_fragment_data = {
        'name': remover.name,
        'attypes': remover.attypes,
        'bond': remover.bond,
        'coord': remover.coord,
        'long_bond': remover.long_bond_idx,
        'free_atom': remover.free_atom
    }
    return new_fragment_data


def synthesize_two_fragment(fragment1_data, fragment2_data):
    """Synthesize fragment1 and fragment2 data, which have one free atoms for each fragment.
    Each Molecule data have
        'name'(str),
        'attypes'(List[str]),
        'bond'(List[List[int]]),
        'long_bond'(List[List[int]]),
        'coord'(List[List[float]]),
        'free_atom'(List[int]).

    :param fragment1_data: Fragment1 Data comprise attypes, bond, coord, and free_atom.
    :param fragment2_data: Fragment2 Data comprise attypes, bond, coord, and free_atom.
    :return: new_mol_data(dict): New molecule data comprise attypes, bond, coord.
    """
    synthesizer = chemical_synthesis.Synthesis(fragment1_data, fragment2_data)
    synthesizer.synthesize_name()
    synthesizer.synthesize_attypes()
    synthesizer.synthesize_coord()
    synthesizer.arrange_atom_index_of_bond()
    synthesizer.arrange_atom_index_of_free_atom()
    synthesizer.synthesize_bond()
    synthesizer.synthesize_free_atom()
    synthesizer.arrange_bond_index_of_long_bond()
    synthesizer.synthesize_long_bond()
    synthesizer.arrange_coord_to_place_fragment()
    synthesizer.add_bond()
    if not synthesizer.free_atom_new_mol:
        new_mol_data = {
            'name': synthesizer.name_new_mol,
            'attypes': synthesizer.attypes_new_mol,
            'bond': synthesizer.bond_new_mol,
            'long_bond': synthesizer.long_bond_idx_new_mol,
            'coord': synthesizer.coord_new_mol,
            'free_atom': synthesizer.free_atom_new_mol
        }
        return new_mol_data
