import sys, os

from molecule_synthesizer.models import fragment

sys.path.append(os.pardir)

import pytest

from molecule_synthesizer.models.fragment import *


@pytest.fixture
def test_fragment() -> Fragment:
    fragment = Fragment('test_fragment', new_molecule=True)
    fragment.attypes = ['', 'C', 'C', 'N', 'H', 'H']
    fragment.bond = [[1, 2], [1, 3], [1, 4], [2, 5]]
    fragment.long_bond_idx = [2, 3]
    fragment.coord= [[], [0., 0., 0.], [1., 1., 1.], [2., 2., 2.], [3., 3., 3.], [4., 4., 4.]]
    fragment.bind_fragment = ['', 'benzothiazole', 'aryl', '', '', '']
    return fragment


class TestRemoveHydrogenTestFragmentBenzothiazole(object):
    @pytest.fixture
    def remove_hydrogen_test_fragment(self, test_fragment) -> None:
        test_fragment.remove_hydrogen('benzothiazole')

    def test_attypes(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.attypes == ['', 'C', 'C', 'N', 'H']

    def test_bond(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.bond == [[1, 2], [1, 3], [2, 4]]

    def test_long_bond_idx(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.long_bond_idx == [2]

    def test_coord(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.coord == [[], [0., 0., 0.], [1., 1., 1.], [2., 2., 2.], [4., 4., 4.]]

    def test_bind_fragment(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.bind_fragment == ['', 'benzothiazole', 'aryl', '', '']

    def test_free_atom(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.free_atom == [1]


class TestRemoveHydrogenTestFragmentAryl(object):
    @pytest.fixture
    def remove_hydrogen_test_fragment(self, test_fragment) -> None:
        test_fragment.remove_hydrogen('aryl')

    def test_attypes(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.attypes == ['', 'C', 'C', 'N', 'H']

    def test_bond(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.bond == [[1, 2], [1, 3], [1, 4]]

    def test_long_bond_idx(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.long_bond_idx == [2]

    def test_coord(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.coord == [[], [0., 0., 0.], [1., 1., 1.], [2., 2., 2.], [3., 3., 3.]]

    def test_bind_fragment(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.bind_fragment == ['', 'benzothiazole', 'aryl', '', '']

    def test_free_atom(self, test_fragment, remove_hydrogen_test_fragment):
        assert test_fragment.free_atom == [2]


class DummyFragment():
    def __init__(self, name: str = None, new_molecule: bool = False):
        self.name = name

        self.attypes = []
        self.long_bond_idx = []
        self.bond = []
        self.coord = []
        self.bind_fragment = []
        self.free_atom = []


class TestSynthesisAllFragments(object):
    @pytest.fixture
    def synthesizer(self, monkeypatch) -> Synthesis:
        fragments = {'benzothiazole': 'test_benzothiazole',
                     'amide': 'test_amide',
                     'aryl': 'test_aryl',
                     'alcohol1': 'test_alcohol1',
                     'alcohol2': 'test_alcohol2',
                     'modifier': 'test_modifier'
                     }
        monkeypatch.setattr(Fragment, '__init__', DummyFragment.__init__)
        synthesizer = fragment.Synthesis(fragments)
        return synthesizer

    @pytest.fixture
    def set_attypes(self, synthesizer, monkeypatch) -> None:
        monkeypatch.setattr(synthesizer.benzothiazole, 'attypes', ['', 'B'])
        monkeypatch.setattr(synthesizer.amide, 'attypes', ['', 'Am'])
        monkeypatch.setattr(synthesizer.aryl, 'attypes', ['', 'Ar'])
        monkeypatch.setattr(synthesizer.alcohol1, 'attypes', ['', 'Al', 'A1'])
        monkeypatch.setattr(synthesizer.alcohol2, 'attypes', ['', 'Al', 'A2'])
        monkeypatch.setattr(synthesizer.modifier, 'attypes', ['', 'Cl'])

    @pytest.fixture
    def set_bond(self, synthesizer, monkeypatch) -> None:
        monkeypatch.setattr(synthesizer.benzothiazole, 'bond', [[1, 2]])
        monkeypatch.setattr(synthesizer.amide, 'bond', [[1, 2]])
        monkeypatch.setattr(synthesizer.aryl, 'bond', [[1, 2]])
        monkeypatch.setattr(synthesizer.alcohol1, 'bond', [[1, 2]])
        monkeypatch.setattr(synthesizer.alcohol2, 'bond', [[1, 2]])
        monkeypatch.setattr(synthesizer.modifier, 'bond', [[1, 2]])

    @pytest.fixture
    def set_bind_fragment(self, synthesizer, monkeypatch) -> None:
        monkeypatch.setattr(synthesizer.benzothiazole, 'bind_fragment', ['', ''])
        monkeypatch.setattr(synthesizer.benzothiazole, 'bind_fragment', ['', ''])
        monkeypatch.setattr(synthesizer.benzothiazole, 'bind_fragment', ['', ''])
        monkeypatch.setattr(synthesizer.benzothiazole, 'bind_fragment', ['', ''])
        monkeypatch.setattr(synthesizer.benzothiazole, 'bind_fragment', ['', ''])

    def test_synthesize_attypes(self, synthesizer, monkeypatch, set_attypes):
        synthesizer.synthesize_attypes()
        assert synthesizer.new_molecule.attypes == ['', 'B', 'Am', 'Ar', 'Al', 'A1', 'Al', 'A2', 'Cl']

    def test_synthesize_bond(self, synthesizer, monkeypatch):

    def test_get_adding_bond(self, synthesizer, set_bind_fragment):
        adding_bond = synthesizer.get_adding_bond(synthesizer.benzothiazole, 4)
        assert adding_bond == []



