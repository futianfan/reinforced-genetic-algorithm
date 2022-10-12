# from tdc import utils
# names = utils.retrieve_benchmark_names('Docking_Group')
# print(names)
# pyscreener_path = '/project/molecular_data/graphnn/pyscreener/'
# from tdc.benchmark_group import docking_group
# group = docking_group(path = 'data/', 
#                 file_format='1iep_docking', 
#                 pyscreener_path = pyscreener_path)

# benchmark = group.get('DRD3', num_max_call = 1000) 

import numpy as np 
import os 
import yaml 
from tdc import Evaluator, Oracle 
div_evaluator = Evaluator(name = 'Diversity')
nov_evaluator = Evaluator(name = 'Novelty')
qed_evaluator = Oracle('qed')
sa_evaluator = Oracle('sa')

base_path = '/project/molecular_data/graphnn/mol_opt/main'
method_list = ['graph_ga', 'jt_vae', 'MARS', 'moldqn', 'rationaleRL', 'REINVENT', 'screening', 'selfies_ga', ]
# method_list = ['smiles_ga']
def _normalize_docking_score(raw_score):
	return 1/(1+np.exp((raw_score+7.5)))

def reverse_normalize(normalize_score): 
	return np.log(1/normalize_score - 1)-7.5

#### test ####
# raw_score_list = [-20, -15, -10, -5, 0, 5] 
# for raw_score in raw_score_list:
# 	print(raw_score, reverse_normalize(_normalize_docking_score(raw_score))) 

def evaluate_from_yaml_file(yaml_file):
	""" 
	Args: 
		- yaml_file

	Return:
		- top-1
		- top-10 
		- top-100 
		- novelty 
		- diversity 
		- qed 
		- sa 
	"""
	result = yaml.load(open(yaml_file, "r").read(), Loader = yaml.Loader)
	result_lst = [(smiles, reverse_normalize(normalize_score), idx) for smiles, (normalize_score, idx) in result.items()]
	top_100 = np.mean([i[1] for i in result_lst[:100]])
	top_10 = np.mean([i[1] for i in result_lst[:10]])
	top_1 = result_lst[0][1]
	top_100_smiles = [i[0] for i in result_lst[:100]]
	div = div_evaluator(top_100_smiles)
	qed_score = np.mean(qed_evaluator(top_100_smiles))
	sa_score = np.mean(sa_evaluator(top_100_smiles))
	return top_100, top_10, top_1, div, qed_score, sa_score 



for method in method_list:
	print('-------- ' + method + ' ----------')
	result_path = os.path.join(base_path, method, 'results')
	files = list(os.listdir(result_path))
	files = list(filter(lambda x:'docking' in x, files))
	result_lst = []
	for file in files:
		file = os.path.join(result_path, file)
		result = evaluate_from_yaml_file(file)
		result_lst.append(result)
	top_100 = np.mean([i[0] for i in result_lst]), np.std([i[0] for i in result_lst]), 
	top_10 = np.mean([i[1] for i in result_lst]), np.std([i[1] for i in result_lst]), 
	top_1 = np.mean([i[2] for i in result_lst]), np.std([i[2] for i in result_lst]), 
	div = np.mean([i[3] for i in result_lst]), np.std([i[3] for i in result_lst]), 
	qed_score = np.mean([i[4] for i in result_lst]), np.std([i[4] for i in result_lst]), 
	sa_score = np.mean([i[5] for i in result_lst]), np.std([i[5] for i in result_lst]), 
	print('\ttop 100', top_100)
	print('\ttop 10', top_10)
	print('\ttop 1', top_1)
	print('\tdiv', div)
	print('\tqed', qed_score)
	print('\tsa', sa_score)






