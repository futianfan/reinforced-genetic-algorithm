from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


with open('mutation_example.txt', 'r') as fin:
	lines = fin.readlines() 

idx = 0
for line in lines[idx:idx+1]:
	input_smiles, smart, output_smiles = line.split()[:3]
	mol = Chem.MolFromSmiles(input_smiles, sanitize=False)
	Draw.MolToFile(mol, 'figure/mutation_input.png', )
	# try:
	# if True: 
	# 	mol = Chem.MolFromSmarts(smart)
	# 	d2d = Draw.MolDraw2DSVG(250,200)
	# 	d2d.DrawMolecule(mol)
		# Draw.MolToFile(mol, 'figure/mutation_smart.png', )
	# except:
		# pass 

	#### draw SMARTS
	# https://www.rdkit.org/docs/GettingStartedInPython.html
	print(smart)
	rxn = AllChem.ReactionFromSmarts(smart)
	d2d = Draw.MolDraw2DCairo(800,300)
	d2d.DrawReaction(rxn)
	png = d2d.GetDrawingText()
	open('./figure/mutation_smart.png','wb+').write(png) 


	mol = Chem.MolFromSmiles(output_smiles, sanitize=False)
	Draw.MolToFile(mol, 'figure/mutation_output.png', )


# >>> from rdkit.Chem import Draw
# >>> rxn = AllChem.ReactionFromSmarts('[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])O.[N-:13]=[N+:14]=[N-:15]>C(Cl)Cl.C(=O)(C(=O)Cl)Cl>[cH:5]1[cH:6][c:7]2[cH:8][n:9][cH:10][cH:11][c:12]2[c:3]([cH:4]1)[C:2](=[O:1])[N:13]=[N+:14]=[N-:15]',useSmiles=True)
# >>> d2d = Draw.MolDraw2DCairo(800,300)
# >>> d2d.DrawReaction(rxn)
# >>> png = d2d.GetDrawingText()
# >>> open('./images/reaction1.o.png','wb+').write(png) 

# import json
# reaction_file = "autogrow/operators/mutation/smiles_click_chem/reaction_libraries/all_rxns/All_Rxns_rxn_library.json"
# smiles_file = 'source_compounds/naphthalene_smiles.smi'

# with open(smiles_file, 'r') as fin:
# 	smiles_lst = fin.readlines() 
# smiles_lst = [line.split()[0] for line in smiles_lst]

# mol_list = []
# for smiles in smiles_lst:
# 	mol = Chem.MolFromSmiles(smiles, sanitize=False)
# 	mol_list.append(mol) 

# reaction_dict = json.load(open(reaction_file))
# # print(reaction_dict)

# for k,v in reaction_dict.items(): 
#     a_reaction_dict = v 
#     rxn = AllChem.ReactionFromSmarts(str(a_reaction_dict["reaction_string"]))
#     rxn.Initialize()
#     for mol in mol_list:
#         try:
#             y = rxn.RunReactants((mol,))
#             print(y)
#         except:
#             pass 

