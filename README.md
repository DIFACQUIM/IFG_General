# IFG_General
Repository for the work Functional group and diversity analysis of BIOFACQUIM: A Mexican natural product database.

Frgament based implementation of

    An algorithm to identify functional groups in organic molecules
    Peter Ertl

    https://jcheminf.springeropen.com/articles/10.1186/s13321-017-0225-z
    
Based on the work described in
    
    https://github.com/rdkit/rdkit/tree/master/Contrib/IFG
    By:
    Richard Hall,
    Guillaume Godin

Usage:
    Python
    # Returns a list of the identified functional groups in a molecule as fragments with environment carbons replaced by dummy atoms (* for aliphatic and *1 for aromatic)
    
    fgs = identify_funtional_groups(smiles)
    print(fgs)
    
Example:
    fgs = identify_funtional_groups('Clc1ccc2Oc3ccccc3N=C(N4CCNCC4)c2c1')
    print(fgs)
    ['[*][Cl]', '[*][O][*]', '[*][N]=[C]([*])[N]([*])[*]', '[*][N]([*])[H]']
        
    
    
