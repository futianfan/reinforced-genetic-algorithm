from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolHash

# smiles0 = "C1CC2=C3C(=CC=C2)C(=CN3C1)[C@H]4[C@@H](C(=O)NC4=O)C5=CNC6=CC=CC=C65"
smiles0 = "C1CC2=C3C(=CC=C2)C(=CN3C1)[C]4[C](C(=O)NC4=O)C5=CNC6=CC=CC=C65"
smiles1 = "C1CC2=C3C(=CC=C2)N(CN3C1)C4C(Cl)(C(=O)NC4=O)"
smiles2 = "C4(S)C(C(=O)NC4=O)C5=CNC6=CC=CC=C65"

smiles3 = "C4(S)C(Cl)(C(=O)NC4=O)"


# ligand_new_smiles = smiles_merge.run_main_smiles_merge(vars, selected_smiles_1, smiles_2)

mol0 = Chem.MolFromSmiles(smiles0)
mol1 = Chem.MolFromSmiles(smiles1)
mol2 = Chem.MolFromSmiles(smiles2)
mol3 = Chem.MolFromSmiles(smiles3)
# im = Chem.Draw.MolToImage(mol)

# Draw.MolsToGridImage([mol])
# Draw.MolToFile(mol, 'figure/example1.png', )
Draw.MolToFile(mol0, 'figure/example0.png', )
Draw.MolToFile(mol1, 'figure/example1.png', )
Draw.MolToFile(mol2, 'figure/example2.png', )
Draw.MolToFile(mol3, 'figure/example3.png', )




# ./autogrow/operators/mutation/smiles_click_chem/reaction_libraries/click_chem_rxns/ClickChem_rxn_library.json




