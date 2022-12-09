

from enum import Enum
from typing import Iterable
from rdkit.Chem.rdchem import RWMol, Atom, BondType
from rdkit.Chem import MolToSmiles
from examples import molecule_examples

# All supported tokens
# TODO: Full support of SMILES syntax 

atoms: tuple[str] = ('C', 'N', 'O', 'F', 'Cl', "Br")
branch: tuple[str] = ('(', ')')
ring: tuple[str] = tuple('123456789')
bond: tuple[str] = ('-', '=')
all_tokens: tuple[str] = atoms + branch + ring + bond

# Map from atom token to atomic number
atom_map = {
    'C' : 6,
    'F' : 9,
    'O' : 8,
    'N' : 7,
    'Br': 35,
    'Cl': 17
}

"""
Wrapper Class around an RWMol instance.
"""
class Molecule():

    def __init__(self):
        self._mol = RWMol()

    def add_atom(self, atom: str) -> None:
        self._mol.AddAtom(Atom(atom_map[atom]))

    def add_bond(self, atom1_idx: int, atom2_idx: int, bond_type: BondType) -> None:
        self._mol.AddBond(atom1_idx, atom2_idx)
        bond = self._mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
        bond.SetBondType(bond_type)

    def to_smiles(self) -> str:
        return MolToSmiles(self._mol)

class SmileType(Enum):
    Atom = 1
    Branch = 2
    Ring = 3
    Bond = 4

"""
Maps valid tokens to a given SmileType

    Raises:
        ValueError on unsupported tokens
"""
def _token_to_type(token: str) -> SmileType:
    if token in atoms:
        return SmileType.Atom
    if token in branch:
        return SmileType.Branch
    if token in ring:
        return SmileType.Ring
    if token in bond:
        return SmileType.Bond

    raise ValueError(f'Token {token} is not a supported value')

"""
Splits a smiles string into valid tokens

    Raises:
        ValueError if an unsupported token is encountered.  
"""
def _tokenise(smiles_string: str) -> Iterable[str]:
    tokens = []
    current_token = ''

    # Reverses the string as tokens are uniquely identifyable in
    # reverse order e.g. C and Cl are ambiguous when parsing from
    # left to right. 

    for char in reversed(smiles_string):
        if char.islower():
            current_token += char
        else:
            current_token = char + current_token

            # Unsupported token
            if current_token not in all_tokens:
                raise ValueError(f'Token {current_token} is not a supported value')

            tokens.append(current_token)
            current_token = ''

    tokens.reverse()
    return tokens

"""
Convert a stream of tokens into a Molecule

    Raises:
        ValueError on ill formed string.   
"""
def _parse(tokens: Iterable[str]) -> Molecule:

    mol = Molecule()

    previous_atom_index = None

    # map from ring number to start atom
    rings = {}
    # stack to track start of branches
    branches = []

    current_atom_index = 0
    default_bond_type=BondType.SINGLE
    current_bond_type=default_bond_type

    for token in tokens:
        token_type = _token_to_type(token)
        match token_type:
            case SmileType.Atom:
                mol.add_atom(token)
                if previous_atom_index is not None:
                    mol.add_bond(previous_atom_index, current_atom_index, current_bond_type)
                    current_bond_type = default_bond_type
                previous_atom_index = current_atom_index
                current_atom_index += 1

            case SmileType.Branch:
                if previous_atom_index is None:
                    raise ValueError("Expected atom before branch")
                match token:
                    case '(':
                        branches.append(previous_atom_index)
                    case ')':
                        previous_atom_index = branches.pop()
                
            case SmileType.Ring:
                if previous_atom_index is None:
                    raise ValueError("Expected atom before ring")              

                if token in rings.keys():
                    ring_start_index = rings[token]
                    mol.add_bond(ring_start_index, previous_atom_index, default_bond_type)
                    rings.pop(token)
                else:
                    rings[token] = previous_atom_index

            case SmileType.Bond:
                if previous_atom_index is None:
                    raise ValueError("Expected atom before bond")
                match token:
                    case '-':
                        current_bond_type = BondType.SINGLE
                    case '=':
                        current_bond_type = BondType.DOUBLE

    # TODO: Check if all rings and branches are closed at end of string
    # TODO: Check final token is sensible (no bond/branch start as last token?)

    return mol 

"""
Create a molecule from a valid smiles string.

Note: Only a subset of smiles syntax is supported.
"""
def smiles_to_molecule(smiles_str: str) -> Molecule:

    if not isinstance(smiles_str, str):
        raise TypeError(f'Smiles string required as input')

    return _parse(_tokenise(smiles_str))
