import sys, os
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
