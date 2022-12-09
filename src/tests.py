from typing import Iterable
from parser import SmileType, _tokenise, smiles_to_molecule, _token_to_type
from examples import molecule_examples
from rdkit.Chem import MolToSmiles, MolFromSmiles
import pytest

# Tests for public interface smiles_to_molecule
# TODO: Add more comprehensive testing

def test_smiles_to_molecule_produces_correct_molecular_representation() -> None:
    # from examples.py
    for example in molecule_examples:
        mol = smiles_to_molecule(example)
        # Use equality with rdkit's smiles parser to verify the output
        assert(mol.to_smiles() == MolToSmiles(MolFromSmiles(example)))

def test_smiles_to_molecule_on_invalid_input_type() -> None:
    with pytest.raises(TypeError) as e:
        smiles_to_molecule(9)
    assert e.match(f'Smiles string required as input')

def test_smiles_to_molecule_on_empty_str() -> None:
    mol = smiles_to_molecule('')
    assert mol.to_smiles() is ''

def test_smiles_to_molecule_on_invalid_smiles_str() -> None:
    with pytest.raises(ValueError) as e:
        mol = smiles_to_molecule('==-')
    assert e.match("Expected atom before bond")

def test_smiles_to_molecule_on_invalid_token_in_str() -> None:
    with pytest.raises(ValueError) as e:
        mol = smiles_to_molecule('==-:')
    assert e.match(f'Token : is not a supported value')

# Test _token_to_type

@pytest.mark.parametrize(
    ("token", "expected_type"),
    [   
        # TODO: Add rest of token types
        ('C', SmileType.Atom),
        ('Cl', SmileType.Atom),
        ('F', SmileType.Atom),
        ('O', SmileType.Atom),
        ('N', SmileType.Atom),
        ('Br', SmileType.Atom),
    
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
    
# Test _tokenise

@pytest.mark.parametrize(
    ("input_str", "expected_tokens"),
    [   
        ('C', ['C']),
        ('CCl', ['C','Cl']),
        ('CCl=C', ['C', 'Cl', '=', 'C'])
    ]
)
def test_tokenise_produces_correct_tokens(input_str: str, expected_tokens: Iterable[str]) -> None:
    assert expected_tokens == _tokenise(input_str)

# TODO: Invalid tests