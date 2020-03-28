# This function is used to standardize compounds in SMILES form and it works as follows: 
    # Compounds not parseable by RDKit returns "Error 1"
    # The largest component of multiple component compounds is retained. 
    # Compounds consisting of any element other than H, B, C, N, O, F, Si, P, S, Cl, Se, Br and I returns "Error 2"
    # Compounds without the mentioned errors are neutralized and reionized to subsequently generate a canonical tautomer.
    # Any other problem returns "Check manually"

from rdkit import Chem
from molvs.standardize import Standardizer
from molvs.charge import Uncharger
from molvs.charge import Reionizer
from molvs.fragment import LargestFragmentChooser
from molvs.tautomer import TautomerCanonicalizer
from rdkit.Chem.rdmolops import RemoveStereochemistry

def FullStandardization(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol == None:
            # If rdkit could not parse the smiles, returns Error 1
            return "Error 1"
        else:            
            STD = Standardizer()
            LFC = LargestFragmentChooser()
            UC = Uncharger()
            RI = Reionizer()
            TC = TautomerCanonicalizer()
            
            mol = STD(mol)
            mol = LFC(mol)
            
            allowed_elements = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"}
            actual_elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
            if len(actual_elements-allowed_elements) == 0:
                mol = UC(mol)
                mol = RI(mol)
                RemoveStereochemistry(mol)
                mol = TC(mol)
                return Chem.MolToSmiles(mol)
            else:
                # If molecule contains other than the allowed elements, returns "Error 2"
                return "Error 2"
    except:
        return "Check manually"   
