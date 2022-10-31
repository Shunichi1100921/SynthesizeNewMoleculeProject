import datetime
import pathlib
import sys, os
from typing import Dict

import pytest

from molecule_synthesizer.models.fragment import Fragment

sys.path.append(os.pardir)

from molecule_synthesizer.controller import synthesis


@pytest.fixture
def get_base_dir_path():
    """
    return: path: new_molecule_DDMMM_YYYY
    """

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    today = datetime.date.today().strftime('%d%b_%Y')
    data_base_dir_path = os.path.join(base_dir, 'Data', f'new_molecule_{today}')

    return data_base_dir_path


@pytest.fixture
def test_full_fragments() -> Dict[str, str]:
    fragments = {'benzothiazole': 'F4',
                 'amide': 'F1',
                 'aryl': 'F5',
                 'alcohol1': 'F2',
                 'alcohol2': 'F25',
                 'modifier': 'F13'}
    return fragments


@pytest.fixture
def test_loss_benzothiazole() -> Dict[str, str]:

    fragments = {'amide': 'F1',
                 'aryl': 'F5',
                 'alcohol1': 'F2',
                 'alcohol2': 'F25',
                 'modifier': 'F13'}
    return fragments


@pytest.fixture
def test_loss_amide() -> Dict[str, str]:

    fragments = {'benzothiazole': 'F4',
                 'aryl': 'F5',
                 'alcohol1': 'F2',
                 'alcohol2': 'F25',
                 'modifier': 'F13'}
    return fragments


@pytest.fixture
def test_loss_aryl() -> Dict[str, str]:

    fragments = {'benzothiazole': 'F4',
                 'amide': 'F1',
                 'alcohol1': 'F2',
                 'alcohol2': 'F25',
                 'modifier': 'F13'}
    return fragments


@pytest.fixture
def test_loss_alcohol1() -> Dict[str, str]:

    fragments = {'benzothiazole': 'F4',
                 'amide': 'F1',
                 'aryl': 'F5',
                 'alcohol2': 'F25',
                 'modifier': 'F13'}
    return fragments


@pytest.fixture
def test_loss_alcohol2() -> Dict[str, str]:

    fragments = {'benzothiazole': 'F4',
                 'amide': 'F1',
                 'aryl': 'F5',
                 'alcohol1': 'F2',
                 'modifier': 'F13'}
    return fragments


@pytest.fixture
def test_loss_modifier() -> Dict[str, str]:

    fragments = {'benzothiazole': 'F4',
                 'amide': 'F1',
                 'aryl': 'F5',
                 'alcohol1': 'F2',
                 'alcohol2': 'F25',}
    return fragments


@pytest.fixture
def test_loss_alcohol1_alcohol2() -> Dict[str, str]:

    fragments = {'benzothiazole': 'F4',
                 'amide': 'F1',
                 'aryl': 'F5',
                 'modifier': 'F13'}
    return fragments


@pytest.fixture
def test_loss_alcohol1_alcohol2_modifier() -> Dict[str, str]:

    fragments = {'benzothiazole': 'F4',
                 'amide': 'F1',
                 'aryl': 'F5'}

    return fragments


def test_synthesizer_for_machineLearning_testMolecule(test_full_fragments, get_base_dir_path):
    synthesis.synthesize_molecule_for_machine_learning(test_full_fragments)
    new_mol_file_path = os.path.join(get_base_dir_path, 'F4_F1_F5_F2_F25_F13')
    assert os.path.exists(new_mol_file_path)


def test_synthesizer_for_machineLearning_loss_benzothiazole(test_loss_benzothiazole):
    with pytest.raises(Exception) as e:
        _ = synthesis.synthesize_molecule_for_machine_learning(test_loss_benzothiazole)

    assert str(e.value) == 'Lack necessary fragment.'

def test_synthesizer_for_machineLearning_loss_amide(test_loss_amide):
    with pytest.raises(Exception) as e:
        _ = synthesis.synthesize_molecule_for_machine_learning(test_loss_amide)

    assert str(e.value) == 'Lack necessary fragment.'


def test_synthesizer_for_machineLearning_loss_aryl(test_loss_aryl):
    with pytest.raises(Exception) as e:
        _ = synthesis.synthesize_molecule_for_machine_learning(test_loss_aryl)

    assert str(e.value) == 'Lack necessary fragment.'


def test_synthesizer_for_machineLearning_loss_alcohol1(test_loss_alcohol1):
    with pytest.raises(Exception) as e:
        _ = synthesis.synthesize_molecule_for_machine_learning(test_loss_alcohol1)

    assert str(e.value) == 'Lack necessary fragment.'


def test_synthesizer_for_machineLearning_loss_alcohol2(test_loss_alcohol2, get_base_dir_path):
    new_mol_file_path = os.path.join(get_base_dir_path, 'F4_F1_F5_F2_F13')
    assert os.path.exists(new_mol_file_path)


def test_synthesizer_for_machineLearning_loss_alcohol1_alcohol2(test_loss_alcohol1_alcohol2, get_base_dir_path):
    new_mol_file_path = os.path.join(get_base_dir_path, 'F4_F1_F5_F13')
    assert os.path.exists(new_mol_file_path)


def test_synthesizer_for_machineLearning_loss_alcohol1_alcohol2_modifier(test_loss_alcohol1_alcohol2_modifier,
                                                                         get_base_dir_path):
    new_mol_file_path = os.path.join(get_base_dir_path, 'F4_F1_F5')
    assert os.path.exists(new_mol_file_path)
