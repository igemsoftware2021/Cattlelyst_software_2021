
#!/usr/bin/env python
"""
Author: Sara Benito . Master Thesis 2018- Metabolic modelling of Actinobacteria
The scripts calls Optknock_robustknock.py that runs the algorithms opt_knock and Robust_knock to find 
possible knock-out that could lead to an over- production of the target products: cinnamaldehyde, benzoic acid and curcumin 
Cobrapy and functions are used from https://github.com/opencobra/m_model_collection/
This requires [cobrapy](https://opencobra.github.io/cobrapy) latest version (14.1, Revision b1aca5bc)
Python 2.7 is used as the programming language
"""

# Import statements

import cobra
import os
import warnings
import re
from itertools import chain
from time import time
import cobra.test
import sympy
import scipy
import scipy.io

import pandas
from cobra import Model, Reaction, Metabolite

import numpy as np
import pandas as pd
from cobra.util.solver import check_solver_status


from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)


import cobra as cb
from copy import deepcopy
from pulp import LpProblem, LpMaximize, LpMinimize, LpVariable, LpAffineExpression, \
                 LpContinuous, LpBinary, LpStatusOptimal, lpSum, LpStatus

#  Optknock_robustknock 
                 
from pipeline.scripts.Optknock_robustknock import OptAndRob
from pipeline.scripts.analysis import get_biomass_equation # added by Delielena Poli in date 20/09/2020 for compatibility with pipeline

def list_excluded_reactions(model):
	"""
	Define a list of reactions that can be knocked-out. It excludes
	exchange reactions, reactions with no genes, essential reactions and reactions
	with essential genes. Objective function is set for biomass production to study if there
	is growth after single deletion strategy
	
	Input: model,  cobrapy model structure
	
	Output: list of reactions to be knocked-out
	"""
   
	essential_reactions=[]
	essential_genes=[]
	possible_reactions=[]
	possible_reactionsygenes=[]
	null_genes=[]
	
	biomass = get_biomass_equation(model) # modified by Delielena Poli 20/09/2020 for compatibility with pipeline
	model.objective= biomass # modified by Delielena Poli 20/09/2020 for compatibility with pipeline

	#model.objective= 'EX_biomass' # modified by Delielena Poli 20/09/2020 for compatibility with pipeline
	#model.reactions.EX_glyc.lower_bound=-1 # modified by Delielena Poli 20/09/2020 for compatibility with pipeline
	#model.reactions.EX_glc.lower_bound=0 # modified by Delielena Poli 20/09/2020 for compatibility with pipeline

	# #model.reactions.get_by_id('EX_o2').lower_bound=0.
	# #smodel.reactions.get_by_id('EX_o2').upper_bound=0.
	
	# Calculating essential reactions
	reaction=single_reaction_deletion(model,model.reactions[0:])
	for x in reaction:
		for i in range(len(reaction[x])):
			if reaction.growth[i]<10E-06:
				essential_reactions.append(model.reactions[i]) 
	gene=single_gene_deletion(model,model.genes[0:]) #Calculate essential genes
	for x in gene:
		for i in range(len(gene[x])):
			if gene.growth[i]<10E-06:
				essential_genes.append(model.genes[i])
	
	#print(essential_reactions,essential_genes)
	for i in range(len(model.reactions)):
		if model.reactions[i] not in essential_reactions:
			possible_reactions.append(model.reactions[i])
	
	for i in range(len(possible_reactions)):
		if possible_reactions[i].genes!=frozenset([]): #checks If there is no associated gene
			if possible_reactions[i].genes not in essential_genes: # Exclude the essential genes
				possible_reactionsygenes.append(possible_reactions[i].id)	
    
	#eliminate transporters and ATP synthase from the candidate reactions
	#   Added by Delielena Poli in date 18/12/2020
	new_poss_rg = []
	n=0
	for i in possible_reactionsygenes:
		reac = model.reactions.get_by_id(i)
		if 'diffusion' not in reac.name and 'ATP synthase' not in reac.name:
			n+=1
			new_poss_rg.append(i)
			#print(reac.id, reac.name)
	print('Number of candidate reactions: ', n)

	return new_poss_rg
	
if __name__ == "__main__":
	
	#Reading the model
	
	model=cobra.io.read_sbml_model('modelo_prueba_cinna.xml')
	possible_reactionsygenes=list_excluded_reactions(model)
	#print(possible_reactionsygenes)
	#Define target reaction
	
	targets_reactions=['Ex_Cinald_e']
	model.reactions.EX_glyc.lower_bound=-1
	model.reactions.EX_glc.lower_bound=0
	#possible_reactionsygenes.append('PHEt2rpp')
	#possible_reactionsygenes.append('TYRt2rpp')
	# Run Optknock/Robustknock increasing the number of possible knock-outs

	j=2  #It defines the number of possible knock-out
	while j<=12:
		Opt_R,Rob_R=OptAndRob(model,targets_reactions,j,possible_reactionsygenes)
		j=j+1
		print(Opt_R,Rob_R)
		
	# if no excluded reactions:
	# for j in range(10):
		# Opt_R,Rob_R=OptAndRob(model,targets_reactions,j)
		# print(Opt_R,Rob_R)		

	
