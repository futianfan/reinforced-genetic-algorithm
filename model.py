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
			log_output = F.log_softmax(output)
			prob_output = F.softmax(output)
			return log_output, prob_output.tolist()  
		else: 
			### smiles string 
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
		log_output = F.log_softmax(output)
		prob_output = F.softmax(output)
		return log_output, prob_output.tolist()  




def atom2int(atom):
	atom_list = ['C', 'N', 'S', 'O', 'H', 'unknown']
	if atom in atom_list:
		return atom_list.index(atom)
	return len(atom_list)-1

def pdbtofeature(pdbfile, centers, pocket_size):
	""" centers=(center_x, center_y, center_z); pocket_size=(size_x, size_y, size_z) """
	with open(pdbfile, 'r') as fin:
		lines = fin.readlines() 
	def featurize(line):
		if line.split()[0]!='ATOM':
			return None 
		if int(line.split()[1])!=float(line.split()[1]):
			return None
		center_x, center_y, center_z = centers 
		size_x, size_y, size_z = pocket_size  
		# xx, yy, zz = float(line.split()[6]), float(line.split()[7]), float(line.split()[8])
		xx = float(line[30:38])
		yy = float(line[38:46])
		zz = float(line[46:54])
		# print('>>>> ', xx, yy, zz)
		if xx < center_x-size_x/2 or xx > center_x+size_x/2:
			return None
		# print('<<<< ', xx, yy, zz)
		if yy < center_y-size_y/2 or yy > center_y+size_y/2:
			return None
		# print('++++ ', xx, yy, zz)
		if zz < center_z-size_z/2 or zz > center_z+size_z/2:
			return None
		# print('----- ', xx, yy, zz) 
		atom_type = line.split()[-1] 
		atom_type = atom2int(atom_type)
		coordinates = torch.FloatTensor([xx, yy, zz]).view(1, -1)
		return atom_type, coordinates 
	lines = list(map(featurize, lines))
	features = list(filter(lambda x:x is not None, lines))
	atom_idx = torch.LongTensor([feature[0] for feature in features]).view(1,-1)  #### (1,N)
	mask = torch.ByteTensor([True for feature in features]).view(1,-1)  ##### (1,N)
	positions = torch.cat([feature[1] for feature in features], dim=0)  ##### (N,3)
	return atom_idx, positions, mask 

receptor_info_list = [
    ('4r6e', './pdb/4r6e.pdb', -70.76, 21.82, 28.33, 15.0, 15.0, 15.0), 
    ('3pbl', './pdb/3pbl.pdb', 9, 22.5, 26, 15, 15, 15), 
    ('1iep', './pdb/1iep.pdb', 15.6138918, 53.38013513, 15.454837, 15, 15, 15),
    ('2rgp', './pdb/2rgp.pdb', 16.29212, 34.870818, 92.0353, 15, 15, 15),
    ('3eml', './pdb/3eml.pdb', -9.06363, -7.1446, 55.86259999, 15, 15, 15),
    ('3ny8', './pdb/3ny8.pdb', 2.2488, 4.68495, 51.39820000000001, 15, 15, 15),
    ('4rlu', './pdb/4rlu.pdb', -0.73599, 22.75547, -31.23689, 15, 15, 15),
    ('4unn', './pdb/4unn.pdb', 5.684346153, 18.1917, -7.3715, 15, 15, 15),
    ('5mo4', './pdb/5mo4.pdb', -44.901, 20.490354, 8.48335, 15, 15, 15),
    ('7l11', './pdb/7l11.pdb', -21.81481, -4.21606, -27.98378, 15, 15, 15), ]

receptor2pdbfeature = dict()

for receptor_info in receptor_info_list:
	name_of_receptor, filename_of_receptor, center_x, center_y, center_z, size_x, size_y, size_z = receptor_info
	# print('------ ' + name_of_receptor + ' --------')
	atom_idx, positions, mask = pdbtofeature(pdbfile=filename_of_receptor, 
											 centers=(center_x, center_y, center_z), 
											 pocket_size=(size_x, size_y, size_z))
	receptor2pdbfeature[name_of_receptor] = (atom_idx, positions, mask)


def pdbqtvina2feature(pdbqt_file):
	with open(pdbqt_file, 'r') as fin:
		lines = fin.readlines() 
	lines = [line.strip() for line in lines]
	line_indx = lines.index("MODEL 2")
	lines = lines[1:line_indx]
	def featurize(line):
		if line.split()[0]!='HETATM':
			return None 
		atom = line.split()[2] #### 'C12'
		atom = [i for i in atom if i > '9'] #### 'C'
		atom = ''.join(atom)
		atom = atom2int(atom)
		coordinates = torch.FloatTensor([float(i) for i in line.split()[6:9]]).view(1, -1)
		return atom, coordinates
		### 2, torch.Tensor([-82.905, 15.268, 40.501])
	lines = list(map(featurize, lines))
	features = list(filter(lambda x:x is not None, lines))
	atom_idx = torch.LongTensor([feature[0] for feature in features]).view(1,-1)  #### (1,N)
	mask = torch.ByteTensor([True for feature in features]).view(1,-1)  #### (1,N)
	positions = torch.cat([feature[1] for feature in features], dim=0)  #### (N,3)
	return atom_idx, positions, mask  

def featurize_receptor_and_ligand(name_of_receptor, pdbqt_file):
	# receptor_atom_idx, receptor_positions, receptor_mask = pdbtofeature(pdbfile, centers, pocket_size)
	receptor_atom_idx, receptor_positions, receptor_mask = receptor2pdbfeature[name_of_receptor]

	ligand_atom_idx, ligand_positions, ligand_mask = pdbqtvina2feature(pdbqt_file) 
	atom_idx = torch.cat([receptor_atom_idx, ligand_atom_idx], dim=1) #### (1,N)
	positions = torch.cat([receptor_positions, ligand_positions], dim=0) 
	positions = torch.unsqueeze(positions, 0)  ##### (1,N,3)
	mask = torch.cat([receptor_mask, ligand_mask], dim=1) ###### (1,N)	
	return atom_idx, positions, mask 


def featurize_receptor_and_ligand_list(name_of_receptor, pdbqt_file_list):
	""" TODO """
	# receptor_atom_idx, receptor_positions, receptor_mask = pdbtofeature(pdbfile, centers, pocket_size)
	receptor_atom_idx, receptor_positions, receptor_mask = receptor2pdbfeature[name_of_receptor]
	feature_list = []
	for pdbqt_file in pdbqt_file_list:
		ligand_atom_idx, ligand_positions, ligand_mask = pdbqtvina2feature(pdbqt_file) 
		atom_idx = torch.cat([receptor_atom_idx, ligand_atom_idx], dim=1) #### (1,N)
		positions = torch.cat([receptor_positions, ligand_positions], dim=0) 
		positions = torch.unsqueeze(positions, 0)  ##### (1,N,3)
		mask = torch.cat([receptor_mask, ligand_mask], dim=1) ###### (1,N)
		feature_list.append((atom_idx, positions, mask))
	return feature_list 


"""https://arxiv.org/pdf/2105.09016.pdf 
TODO
	- test equivariance 
"""
class ENN(nn.Module):
	"""Args:
			1. amino's categories
			2. amino's position
	"""
	def __init__(self, latent_dim = 50, device = torch.device('cpu'), is_one_hot = True, layer = 1, vocab_size=6, coordinate_dim = 3):
		super(ENN, self).__init__() 
		self.latent_dim = latent_dim 
		self.layer = layer 
		self.is_one_hot = is_one_hot
		self.vocab_size = vocab_size 
		self.coordinate_dim = 3 
		self.device = device
		self.aggregate = torch.mean  
		if is_one_hot:
			self.node_embedding = nn.Embedding(vocab_size, latent_dim).to(device)
		
		self.phi_e = nn.Sequential(
								nn.Linear(2*self.latent_dim+1, self.latent_dim), 
								nn.Tanh(), 
								nn.Linear(self.latent_dim, self.latent_dim), 
								nn.Tanh()).to(device) ## 2d+1 -> d

		self.phi_x = nn.Sequential(
								nn.Linear(self.latent_dim, self.latent_dim), 
								nn.Tanh(), 
								nn.Linear(self.latent_dim, 1), 
								nn.Tanh(), 
								).to(device)  ### d->1

		self.phi_h = nn.Sequential(
								nn.Linear(self.latent_dim, self.latent_dim), 
								nn.Tanh(), 
								nn.Linear(self.latent_dim, self.latent_dim), 
								nn.Tanh(),
								).to(device) ### d->d

		self.phi_inf = nn.Sequential(
								nn.Linear(self.latent_dim, 1), 
								nn.Sigmoid(),
								).to(device)  ## d->1

		self.output_mlp = nn.Linear(self.latent_dim, 1)

		# self.output_mlp = nn.Sequential(
		# 						nn.Linear(self.latent_dim, 1), 
		# 						nn.Sigmoid(),
		# 						)

	def forward(self, input_data, coordinate, mask):
		"""
		Args:
			input_data: LongTensor(b,N) & FloatTensor(b,N,d)
			coordinate: b,N,3
			mask: b,N
			where b = batchsize, N = max_num_of_atom  
		
		Returns:
			(b,1)

		"""
		transform = False
		H = self.node_embedding(input_data) if self.is_one_hot else input_data
		if H.dim() == 4:  ##### (a1,a2,N,d)
			transform = True 
			a1,a2,a3,a4 = H.shape
			H = H.view(-1,a3,a4) ## b,N,d 
			b1,b2,b3,b4 = coordinate.shape 
			coordinate = coordinate.view(-1,b3,b4) ### b,N,3
			d1,d2,d3 = mask.shape 
			mask = mask.view(-1,d3)  ## b,N

		b, N = H.shape[0], H.shape[1]
		X = coordinate ### b,N,3 
		mask_expand = mask.unsqueeze(-1) #### b,N,1 
		mask_expand2 = mask_expand.permute(0,2,1) ### b,1,N
		mask_square = mask_expand * mask_expand2 ### b,N,N 
		mask_square = mask_square.unsqueeze(-1) ### b,N,N,1 

		for l in range(self.layer):
			### 1. m_ij = phi_e(h_i, h_j, ||x_i^l - x_j^l||^2)
			H1 = H.unsqueeze(2).repeat(1,1,N,1) ### b,N,N,d
			H2 = H.unsqueeze(1).repeat(1,N,1,1) ### b,N,N,d
			x1 = X.unsqueeze(2).repeat(1,1,N,1) ### b,N,N,3
			x2 = X.unsqueeze(1).repeat(1,N,1,1) ### b,N,N,3
			x12 = torch.sum((x1-x2)**2 * mask_square, dim=-1, keepdim=True) ### b,N,N,1
			H12x = torch.cat([H1,H2,x12], -1) ### b,N,N,2d+1
			M = self.phi_e(H12x)*mask_square ### b,N,N,d 
			### 2. e_ij = phi_inf(m_ij)
			E = self.phi_inf(M) ### b,N,N,1
			### 3. m_i = \sum e_ij m_ij
			M2 = torch.sum(M*E,1) ## b,N,d 
			### 4. x_i^{l+1} = x_i^l + \sum_{j\neq i} (x_i^l - x_j^l) phi_x(m_ij) 
			X = X + torch.sum((x1 - x2) * mask_square * self.phi_x(M), dim=1) ## b,N,3
			### 5. h_i^{l+1} = phi_h(h_i^l, m_i)
			H = self.phi_h(M2) + H   ### b,N,d 
			H = H * mask_expand  ### b,N,d  

		if transform:
			H = H.view(a1,a2,a3,a4)
			mask = mask.view(d1,d2,d3)
		H = self.aggregate(H*mask.unsqueeze(-1), dim = -2) 
		H = nn.ReLU()(H)
		H = self.output_mlp(H)
		return H 

	def forward_ligand_list(self, name_of_receptor, pdbqtvina_list):
		feature_list = featurize_receptor_and_ligand_list(name_of_receptor, pdbqtvina_list)
		output_list = []
		for atom_idx, positions, mask in feature_list: 
			output = self.forward(atom_idx, positions, mask) #### [1,1]
			output_list.append(output)
		outputs = torch.cat(output_list, dim=0).view(-1)
		log_output = F.log_softmax(outputs, 0)
		prob_output = F.softmax(outputs, 0)
		# print("output probability", prob_output.tolist()[:3])
		return log_output, prob_output.tolist()  




if __name__ == "__main__":


	# model = Ligand2D() 
	# smiles = ['CCC', 'CCC']
	# output = model(smiles)
	# print(output.shape, output)
	# output = model(smiles[0])

	# model = Ligand2D_product() 
	# output = model(smiles[0], smiles)
	# print(output)

	# atom_idx, positions, mask = pdbqtvina2feature(pdbqt_file='4r6e_example.pdbqt.vina')
	# print(atom_idx, positions, atom_idx.shape, positions.shape)

	# atom_idx, positions, mask = pdbtofeature(pdbfile='./pdb/4r6e.pdb', centers=(-70.76, 21.82, 28.33), pocket_size=(18.0, 18.0, 18.0))

	# print(atom_idx.shape, positions.shape)

	# atom_idx, positions, mask = featurize_receptor_and_ligand(pdbfile='./pdb/4r6e.pdb', 
	# 														  centers=(-70.76, 21.82, 28.33), 
	# 														  pocket_size=(18.0, 18.0, 18.0), 
	# 														  pdbqt_file='4r6e_example.pdbqt.vina')

	# enn = ENN()

	# output = enn(input_data = atom_idx, coordinate = positions, mask = mask)
	# print(output.shape, output)

	# output = enn(input_data = atom_idx, coordinate = positions+2, mask = mask)
	# print(output.shape, output)
	pass 






