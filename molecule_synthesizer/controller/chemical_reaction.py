"""Data: atom types, bond, coord, long-bond, free-atom"""
from molecule_synthesizer.models import chemical_synthesis
from molecule_synthesizer.models.chemical_synthesis import Fragment


def synthesize_two_fragment(fragment1: Fragment, fragment2: Fragment) -> Fragment:
    """Synthesize fragment1 and fragment2, which are Fragment object and have one free atoms for each fragment.
    :param fragment1: Fragment1 Data comprise attypes, bond, coord, and free_atom.
    :param fragment2: Same as fragment1.
    :return: new_mol_data(dict): New molecule data comprise attypes, bond, coord.
    """
    synthesizer = chemical_synthesis.Synthesis(fragment1, fragment2)
    synthesizer.synthesize_attypes()
    synthesizer.synthesize_bond()
    synthesizer.synthesize_free_atom()
    synthesizer.synthesize_long_bond()
    synthesizer.synthesize_coord()
    synthesizer.add_bond()
    return synthesizer.new_molecule
