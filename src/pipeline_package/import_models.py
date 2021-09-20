#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Loading the reference and universal cobra models

This module reads the input file in csv format to get the information
needed for loading the user defined reference and universal BiGG models
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"


import warnings
import os
from os.path import join

import cobra
import cobra.test
from cobra import Model
import csv

warnings.filterwarnings('ignore')

def get_universal_suffix(infile):
    """
    Read input for the type of universal the users wants to use.

    in the input file the suer can select among Gram-positive, 
    Gram-negative, cyanobacteria and Archaea. This can be done by 
    properly filling the universal_reaction_DB column with either
    bacteria, grampos, gramneg, cyanobacteria or archaea respectively. 
    --------------------------------------------------------
    Argument:
    infile--str csv file dictionary like

    Return
    uni: str specific suffix of the name for the universal model
    """
    with open(infile, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            for row in reader:
                if row['universal_reaction_DB'] != '':
                    uni = row['universal_reaction_DB'] 
    return uni

def get_universal_main(data_dir, infile):
    """
    Open universal file with adequate suffix

    Use the suffix obtaiened with get_universal_suffix finction
    to open the universal file form the spacified location. 
    Additionally it modifies the universal model by removing the
    exchange reactions. 
    -------------------------------------------------------
    Arguments:
    data_dir--str path to be read by os.path.join 
    infile--str csv file dictionary like

    Return:
    universal: cobra.Model universal model 
    """
    suffix = get_universal_suffix(infile)
    uni_name = 'universe_'+suffix+'.xml.gz'
    universal = cobra.Model("reaction_reference")
    universal = cobra.io.read_sbml_model(join(data_dir, uni_name)) 
    #TODO: check path of universal 
    
    for react in list(universal.reactions):
        if "EX_" in react.id:
            universal.remove_reactions(react)

    return universal

def get_ID_reference_model(infile):
    """
    Get BiGG ID of the reference model

    The BiGG ID of the reference model must be provided in the input 
    file in correspondance to 'reference_model_ID' column. 
    ----------------------------------------------
    Argument:
    infile--str input file in .csv format dictionary like

    Return
    ref_model: str with BiGG ID of the reference model of a specific 
        organism
    """
    with open(infile, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            for row in reader:
                if row['reference_model_ID'] != '':
                    ref_model = row['reference_model_ID'] 
    
    return ref_model

def get_expression_host(infile):
    """
    Open reference BiGG model

    The BiGG ID of the reference model must be provided in the input 
    file in correspondance to 'reference_model_ID' column. 
    The file has to be already downloaded and saved in the repository
    were the pipeline scripts are.
    ----------------------------------------------
    Argument:
    infile--str input file in .csv format dictionary like

    Return
    exp_host: str expression host species name, correspondant
        to reference model ID. 
    """
    with open(infile, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            for row in reader:
                if row['expression_host'] != '':
                    exp_host = row['expression_host'] 
    
    return exp_host

def get_reference_model(data_dir, infile):
    """
    Open reference BiGG model

    The BiGG ID of the reference model must be provided in the input 
    file in correspondance to 'reference_model_ID' column. 
    The file has to be already downloaded and saved in the repository
    were the pipeline scripts are.
    ----------------------------------------------
    Arguments:
    data_dir--str path to be read by os.path.join
    infile--str input file in .csv format dictionary like

    Return
    model: cobra.Model reference model of a specific organism
    """
    ref_model = get_ID_reference_model(infile)
    
    model = cobra.Model("reference_model")
    model = cobra.io.read_sbml_model(join(data_dir, ref_model))
    return model



    
