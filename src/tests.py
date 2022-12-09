from typing import Iterable
from parser import SmileType, smiles_to_molecule, _token_to_type
from examples import molecule_examples
from rdkit.Chem import MolToSmiles, MolFromSmiles
import pytest

def test_smiles_to_molecule_produces_correct_molecular_representation() -> None:
    # from examples.py
    for example in molecule_examples:
        mol = smiles_to_molecule(example)
        # Use equality with rdkit's smiles parser to verify our output
        assert(mol.to_smiles() == MolToSmiles(MolFromSmiles(example)))

def test_smiles_to_molecule_on_invalid_input() -> None:
    with pytest.raises(TypeError) as e:
        smiles_to_molecule(9)
    assert e.match(f'Smiles string required as input')

# Test _token_to_type

@pytest.mark.parametrize(
    ("token", "expected_type"),
    [
        ('C', SmileType.Atom),
    ]
)
def test_token_to_types_produces_expected_type(token: str, expected_type: SmileType) -> None:
    assert _token_to_type(token) == expected_type

@pytest.mark.parametrize(
    "token",
    [
        'Z', 'Foo', '_', ''
    ]
)
def test_invalid_token_is_not_supported(token: str) -> None:
    with pytest.raises(ValueError, ) as e:
        _token_to_type(token)
    assert e.match(f'Token {token} is not a supported value')
    
#