# Based on the work of Richard Hall and Guillaume Godin
# Additions or modifications are indicated with ##

from rdkit import Chem
from rdkit.Chem import rdmolops

def merge(mol, marked, aset):
    bset = set()
    for idx in aset:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            jdx = nbr.GetIdx()
            if jdx in marked:
                marked.remove(jdx)
                bset.add(jdx)
    if not bset:
        return
    merge(mol, marked, bset)
    aset.update(bset)

## Heteroatoms connected to an aliphatic atom (by single, double or triple bond), to avoid single aromatic heteroatoms
PATT_HETERO = Chem.MolFromSmarts("[!#6;!#1]")
# atoms connected by non-aromatic double or triple bond to any heteroatom
# c=O should not match (see fig1, box 15).  I think using A instead of * should sort that out?
## Using C explicitly indicates non-aromatic carbons, we changed this beacuse we use explicit hydrogens in the molecule, that would be matched by A 
PATT_DOUBLE_TRIPLE = Chem.MolFromSmarts('C=,#[!#6]')
# atoms in non aromatic carbon-carbon double or triple bonds
PATT_CC_DOUBLE_TRIPLE = Chem.MolFromSmarts('C=,#C')
# acetal carbons, i.e. sp3 carbons connected to tow or more oxygens, nitrogens or sulfurs; these O, N or S atoms must have only single bonds
PATT_ACETAL = Chem.MolFromSmarts('[CX4](-[O,N,S])-[O,N,S]')
# all atoms in oxirane, aziridine and thiirane rings
PATT_OXIRANE_ETC = Chem.MolFromSmarts('[O,N,S]1CC1')

PATT_TUPLE = (PATT_HETERO, PATT_DOUBLE_TRIPLE, PATT_CC_DOUBLE_TRIPLE, PATT_ACETAL, PATT_OXIRANE_ETC)

def identify_functional_groups(smi):
    ## We decided to start from a SMILES and add explicit hydrogens inside the function
    mol = Chem.MolFromSmiles(smi)
    mol = rdmolops.AddHs(mol)
    try:
        marked = set()
    ## Since heteroatoms are included in PATT_TUPLE, we remove the first part of the original function 
        for patt in PATT_TUPLE:
            for path in mol.GetSubstructMatches(patt):
                for atomindex in path:
                    marked.add(atomindex)

    #merge all connected marked atoms to a single FG
        groups = []
        while marked:
            grp = set([marked.pop()])
            merge(mol, marked, grp)
            groups.append(grp)
        groups = [list(x) for x in groups]
        
    ## It seems that the initial filtering of heteroatoms was not enough, so we add this to remove groups with only aromatic atoms
        for g in groups:
            group_aromaticity = set([mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in g])
            if group_aromaticity == {True}:
                groups.remove(g)
        
    ## Identify bonds to break and hydrogens to keep for every FG
        bonds = []
        labels = []
        for g in groups:
            group_bonds = []
            group_labels = []
            for idx in g:
                atom = mol.GetAtomWithIdx(idx)
                
                ## Carbon atoms
                if atom.GetAtomicNum() == 6:
                    for nbr in atom.GetNeighbors():
                        ## Carbonyl groups to disciminate between aldehydes and ketones
                        if nbr.GetAtomicNum() == 8 and str(mol.GetBondBetweenAtoms(idx,nbr.GetIdx()).GetBondType()) == "DOUBLE":
                            PreserveH = True
                            break
                        else:
                            PreserveH = False
                    if PreserveH == True:
                        for nbr in atom.GetNeighbors():
                            jdx = nbr.GetIdx()
                            if jdx not in g and nbr.GetAtomicNum() != 1:
                                group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                group_labels.append((0,0))
                    else:
                        for nbr in atom.GetNeighbors():
                            jdx = nbr.GetIdx()
                            if jdx not in g:
                                group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                group_labels.append((0,0))
                ## Nitrogen atoms
                elif atom.GetAtomicNum() == 7:
                    ## To discriminate between anilines and amines (primary, secondary, etc)
                    if len(g) == 1:
                        neigh_atn = [x.GetAtomicNum() for x in atom.GetNeighbors() if x.GetAtomicNum() != 1]
                        if neigh_atn.count(6) == 1:
                            for nbr in atom.GetNeighbors():
                                jdx = nbr.GetIdx()
                                if jdx not in g and nbr.GetAtomicNum() != 1:
                                    group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                    if nbr.GetIsAromatic() == True:
                                        group_labels.append((1,1))
                                    else:
                                        group_labels.append((0,0))
                        else:
                            for nbr in atom.GetNeighbors():
                                jdx = nbr.GetIdx()
                                if jdx not in g and nbr.GetAtomicNum() != 1:
                                    group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                    group_labels.append((0,0))
                    else:
                        for nbr in atom.GetNeighbors():
                            jdx = nbr.GetIdx()
                            if jdx not in g:
                                group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                group_labels.append((0,0))

                ## Oxygen atoms
                elif atom.GetAtomicNum() == 8:
                    ## To discriminate between alcohols from phenols and esthers from carboxylic acids
                    if len(g) == 1:
                        neigh_atn = [x.GetAtomicNum() for x in atom.GetNeighbors() if x.GetAtomicNum() != 1]
                        if len(neigh_atn) == 1 and neigh_atn.count(6) == 1:
                            for nbr in atom.GetNeighbors():
                                jdx = nbr.GetIdx()
                                if jdx not in g and (nbr.GetAtomicNum() != 1):
                                    group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                    if nbr.GetIsAromatic() == True:
                                        group_labels.append((1,1))
                                    else:
                                        group_labels.append((0,0))
                        else:
                            for nbr in atom.GetNeighbors():
                                jdx = nbr.GetIdx()
                                if jdx not in g and nbr.GetAtomicNum() != 1:
                                    group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                    group_labels.append((0,0))                        
                    else:
                        for nbr in atom.GetNeighbors():
                            jdx = nbr.GetIdx()
                            if jdx not in g and nbr.GetAtomicNum() != 1:
                                group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                group_labels.append((0,0))

                ## Sulfur atoms
                elif atom.GetAtomicNum() == 16:
                    if len(g) == 1:
                        for nbr in atom.GetNeighbors():
                            jdx = nbr.GetIdx()
                            if jdx not in g and nbr.GetAtomicNum() != 1:
                                group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                group_labels.append((0,0))
                    else:
                        for nbr in atom.GetNeighbors():
                            jdx = nbr.GetIdx()
                            if jdx not in g:
                                group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                                group_labels.append((0,0))

                else:               
                    for nbr in atom.GetNeighbors():
                        jdx = nbr.GetIdx()
                        if jdx not in g:
                            group_bonds.append(mol.GetBondBetweenAtoms(idx,jdx).GetIdx())
                            group_labels.append((0,0))
            labels.append(group_labels)
            bonds.append(group_bonds)

    ## Build final fragments
        FGS_ENVS = []
        for i in range(len(groups)):
            Frag = Chem.FragmentOnBonds(mol,bonds[i], dummyLabels = labels[i])
            Frags = rdmolops.GetMolFrags(Frag)
            for j in Frags:
                if groups[i][0] in j:
                    FGS_ENVS.append(Chem.MolFragmentToSmiles(Frag, j, canonical=True, allHsExplicit=True))
        FGS_ENVS = list(set(FGS_ENVS))
        for i in FGS_ENVS:
            if Chem.MolFromSmiles(i) == None:
                FG = Chem.MolFromSmarts(i)
            else:
                FG = Chem.MolFromSmiles(i)
            if set([atom.GetIsAromatic() for atom in FG.GetAtoms() if atom.GetSymbol() not in ["*","H"]]) == {True}:
                FGS_ENVS.remove(i)
        return FGS_ENVS 
    
    except:
        ## When the molecules is as small as a single FG
        FGS_ENVS = [Chem.MolToSmiles(mol, canonical=True, allHsExplicit=True)]
        return FGS_ENVS
