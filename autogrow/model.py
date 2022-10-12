import numpy as np 
import torch 
from torch import nn 
import torch.nn.functional as F

from rdkit import Chem, DataStructs 
from rdkit.Chem import AllChem, Descriptors 

def smiles2fp(smiles_string):
	mol = Chem.MolFromSmiles(smiles_string)
	Chem.SanitizeMol(mol)
	fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
	features = np.zeros((1,))
	DataStructs.ConvertToNumpyArray(fp, features)
	fingerprint = torch.from_numpy(features).float().view(1,-1) 	
	return fingerprint ### [1,2048] torch.Tensor

class Ligand2D(nn.Module):
	"""
	input: SMILES
	output: scalar 
	"""
	def __init__(self, ):
		super(Ligand2D, self).__init__()
		self.input_mlp = nn.Linear(2048, 100)
		self.output_mlp = nn.Linear(100, 1)

	def forward(self, smiles_):
		"""
			:param smiles_
				- list of SMILES string
				- SMILES string  
		"""

		if type(smiles_) == list:
			fps = [smiles2fp(s) for s in smiles_]
			fps = torch.cat(fps, 0)
			hidden_state = F.relu(self.input_mlp(fps))
			output = self.output_mlp(hidden_state)
			output = output.view(-1)
			output = F.softmax(output)
			return output 
		else:
			fingerprint = smiles2fp(smiles_)
			hidden_state = F.relu(self.input_mlp(fingerprint))
			output = self.output_mlp(hidden_state)
			return output ### [1,1]


class Ligand2D_product(nn.Module):
	'''
		input:	ligand2d & product_smiles
		output: scalar
	'''
	def __init__(self, ):
		super(Ligand2D_product, self).__init__()
		self.ligand_mlp = nn.Linear(2048, 100)
		self.product_mlp = nn.Linear(2048, 100)
		self.output_mlp = nn.Linear(200, 1)

	def forward(self, ligand_smiles, product_smiles_list):
		n = len(product_smiles_list)
		ligand_fp = smiles2fp(ligand_smiles)
		ligand_embedding = F.relu(self.ligand_mlp(ligand_fp))
		ligand_embedding = ligand_embedding.repeat(n,1)

		product_fps = [smiles2fp(smiles) for smiles in product_smiles_list]
		product_fps = torch.cat(product_fps, 0)
		product_embeddings = F.relu(self.product_mlp(product_fps))

		latent_variable = torch.cat([ligand_embedding, product_embeddings], 1)
		output = self.output_mlp(latent_variable).view(-1)
		output = F.softmax(output)
		return output 





if __name__ == "__main__":
	model = Ligand2D() 
	smiles = ['CCC', 'CCC']
	output = model(smiles)
	print(output.shape, output)
	output = model(smiles[0])

	model = Ligand2D_product() 
	output = model(smiles[0], smiles)
	print(output)




