#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for pathway finding and thermodynamic driving force calculation 

Contains function to define the reactions involved in a pathway by 
using pFBA with the target reaction as objective. The minimal set
of reaction that is identified this way is used for the calculation
of max-min driving force (MDF) values. 
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"


import csv
import numpy as np
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import cobra
from cobra.flux_analysis import phenotype_phase_plane
from cobra.flux_analysis import production_envelope
from cobra.flux_analysis import gapfill, pfba 
import sbtab
from sbtab import SBtab, validatorSBtab, sbml2sbtab
from equilibrator_api import ComponentContribution, Q_
from equilibrator_pathway.pathway import Pathway
from equilibrator_cache.exceptions import (
    MissingDissociationConstantsException, ParseException)
from cobra.exceptions import (
    OPTLANG_TO_EXCEPTIONS_DICT, Infeasible, OptimizationError, SolverNotFound)

from pipeline.scripts.input_parser import *

def get_ex_c_source_metab(input_file, model):
    """
    Identifies the carbon source to use in the optimization problem. 
    
    In the input model there should be only one carbon-source metabolite.
    It should be indicated by adding the BiGG ID to the column 
    --------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace

    Return:
    ex: exchange reaction of the c-source
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        for row in reader:
            if (row['carbon_source'] in ('yes', 'Yes', 'Y', 'y', 'YES') and
                    row['metabolite_BiGG_ID']!= ''):
                metab = row['metabolite_BiGG_ID']
                ex = model.reactions.get_by_id('EX_'+metab+'_e')
    return ex

def get_c_source(input_file):
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            c_source = ''
            for row in reader:
                if row['carbon_source'] == 'yes':
                    c_source += row['metabolite_BiGG_ID']
    return c_source

def get_biomass_equation(model):
    """
    Find the reaction corresponding to the biomass equation

    First retrieves the objective expression which in most models is
    the biomass equation by default. In case the objective has been 
    changed it looks of the reactions having biomass in the name or
    ID.
    --------------------------------------------------------
    Argument:
    model--cobra.Model reference model in BiGG namespace

    Return:
    biomass_r_list[0] or reaction: cobra.Reaction for biomass equation
    """
    expression = str(model.objective.expression)
    y = expression.split()
    z = y[0].split('*')
    id_biomass_equation = z[1]
    biomass_eq = model.reactions.get_by_id(id_biomass_equation)
    if 'iomass' in biomass_eq.id or 'BIOMASS' in biomass_eq.id:
        return biomass_eq
    else:
        biomass_r_list = []
        for i in model.reactions:
            if ('iomass' in i.id or 
                    'BIOMASS' in i.id or 
                    'iomass' in i.name or 
                    'BIOMASS' in i.name):
                biomass_r_list.append(i)
        if len(biomass_r_list)==1:
            return biomass_r_list[0]
        elif (len(biomass_r_list)>1 and 
                'core' in [reaction for reaction in biomass_r_list]):
            return reaction
        else:
            return biomass_r_list[0]

def defining_metabolism(input_file, model):
    """
    Read input file to get type of metabolism

    The field aerobic should be filled with either yes or no
    (or similar strings) to indicate whether oxygen presence
    should be considered in the analysis.
    --------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace

    Return:
    metabolism: string of either aerobic or anaerobic
     in ('yes', 'Yes', 'Y', 'y', 'YES') 
    """
    metabolism = ''
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            for row in reader:
                if row['aerobic']  in ('yes', 'Yes', 'Y', 'y', 'YES'):
                    o2 = model.reactions.get_by_id('EX_o2_e')
                    o2.bounds = (-1000, 1000)
                    metabolism += 'aerobic'
                elif row['aerobic'] in ('No', 'no', 'NO', 'N', 'n'):
                    o2 = model.reactions.get_by_id('EX_o2_e')
                    o2.bounds = (0, 1000)
                    metabolism += 'anaerobic'
                else:
                    continue
    return metabolism


def oxygen_sensitivity_check(model):
    """
    Eliminate oxygen sensitive reactions from the model

    Eliminate reactions which are associated to known 
    oxygen sensitive reactions. 
    ------------------------------------------
    Argument:
    model--cobra.Model reference model in BiGG namespace
    """
    # TODO: remove fdxr containing reactions
    for i in model.reactions:
        #flavodoxin containing reactions
        if 'flxs' in i.reaction:
            #print(i.id, i.name, i.reaction)
            model.remove_reactions([i.id])
        #pyruvate formate lyase
        elif 'PFL' == i.id:
            i.bounds = (0,0)
        #pyruvate synthase
        elif 'POR5' == i.id:
            i.bounds = (0,0)
        #thioredoxin containing reactions
        elif 'trdr' in i.reaction:
            #print(i.id, i.name, i.reaction)
            model.remove_reactions([i.id]) 
        else:
            continue


def set_constraints_c_source(input_file, model):
    """
    Set bounds range of the carbon source 

    Those derive from the coefficient in the balanced equantion. The 
    user must write the coefficient under the column
    coefficient_c_source. 

    The range is obtaiend multiplying the coefficient per 1.1 and 0.9
    to avoid constraining the model too much.  
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            for row in reader:
                if row['coefficient_c_source'] != '':
                    c_reaction = get_ex_c_source_metab(input_file, model)
                    coeff = float(row['coefficient_c_source'])
                    lb_coeff = 1.1 * coeff
                    ub_coeff = 0.9 * coeff
                    c_reaction.bounds = (lb_coeff, ub_coeff)
                else:
                    continue
    

def set_constraints_targets(input_file, model):
    """
    Set bounds range of the target

    Those derive from the coefficient in the balanced equantion. The 
    user must write the coefficient under the column
    coefficient_product.  
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            for row in reader:
                if row['coefficient_product'] != '':
                    #exchange_targets = get_ex_target(input_file, model)
                    #coeff_p = row['coefficient_product'].split()
                    #print(coeff_p)
                    #for i in exchange_targets:
                        #index = exchange_targets.index(i)
                        #i.lower_bound = float(coeff_p[index])
                    coeff_product = row['coefficient_product']
                    exchange_target = get_ex_target(input_file, model)
                    exchange_target.lower_bound = float(coeff_product)
                else:
                    continue


def get_ex_target(input_file, model):
    """
    Get the exchange reaction of the target compound

    Reads the input file to get the ID of the target compound and
    uses it to get the exchange reaction from the reference model.
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace

    Retunr:
    ex: cobra.Reaction, exchange reaction of the target compound
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        #exchanges_targets = []
        for row in reader:
            if row['production'] == 'yes':
                metab = row['metabolite_BiGG_ID']
                ex = model.reactions.get_by_id('EX_'+metab+'_e')
                #exchanges_targets.append(ex)
    return ex


def add_free_balancing_reactions(input_file, model):
    """
    Add reactions for unlimited ATP and cofactors

    To get the minimal set of reactions for the production pathway 
    with pFBA some reaction has to be added to allow cofactor 
    regenaration and unlimited ATP supply.

    Under anaerobic conditions some more cofactors could be used,
    hence if the metabolism is anaerobic additional reactions are 
    included.
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    """
    metabolism_type = defining_metabolism(input_file, model)
    aerobic_free_balancing = [
        'adp_c + h_c + pi_c --> atp_c + h2o_c', 
        'nadp_c + h_c --> nadph_c', 
        'nad_c + h_c --> nadh_c', 
        'fad_c + h_c --> fadh2_c',
        'nadh_c --> nad_c + h_c',
        'q8_c + h_c --> qh2_c',
        'accoa_c --> ac_c + coa_c']
    anaerobic_free_balancing = aerobic_free_balancing + ['fdxox_c + h_c --> fdxrd_c', 
                                'fdxo_2_2_c + h_c --> fdxrd_c', 
                                'fdxo_42_c + h_c --> fdxr_42_c', 
                                'flxso_c + h_c --> flxr_c']
                                #TODO: add thioredoxin ?
    if metabolism_type == 'aerobic':
        for n in range(len(aerobic_free_balancing)):
            if 'endless_'+str(n) not in model.reactions:     
                r = cobra.Reaction('endless_'+str(n))
                model.add_reaction(r)
                r_string = aerobic_free_balancing[n]
                r.build_reaction_from_string(r_string, fwd_arrow='-->')

    elif metabolism_type == 'anaerobic':
        for n in range(len(anaerobic_free_balancing)):
            if 'endless_'+str(n) not in model.reactions:
                r = cobra.Reaction('endless_'+str(n))
                model.add_reaction(r)
                r_string = aerobic_free_balancing[n]
                r.build_reaction_from_string(r_string, fwd_arrow='-->')


def get_raw_reaction_list(input_file, model):
    """
    Run pFBA and get raw minimal set of reactions for production.

    It does pFBA and gets the solutions with non-zero flux. It also 
    excludes the reactions with flux rangeing between 
    -1e-14 and 1e-14.

    Exception:
    if pFBA is infeasible the message is reported and the coefficients
    of substrate and target should be checked. Alternatively, the 
    default boundary reaction of the target (exchange) might not be a
    suitable objective for production, hence other boudary reaction 
    types might be used (demand for instance).
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace

    Return:
    non_zero:list of cobra.Reaction objects with a non zero flux
    """
    biomass = get_biomass_equation(model)
    non_zero = []
    target = get_ex_target(input_file, model)
    substrate = get_ex_c_source_metab(input_file, model)
    model.objective = target
    print('Target: ', target)
    fba_production = model.optimize()
    print('FBA objective value: ', fba_production.objective_value, 
        '\nSubstrate consumption flux: ', fba_production.fluxes[substrate.id], 
        '\nTarget production flux: ', fba_production.fluxes[target.id], 
        '\nBiomass: ', fba_production.fluxes[biomass.id], '\n')
    try:
        pfba_solution = pfba(model)
        print('pFBA objective value: ', pfba_solution.objective_value, 
            '\nSubstrate consumption flux: ', pfba_solution.fluxes[substrate.id], 
            '\nTarget production flux: ', pfba_solution.fluxes[target.id], 
            '\nBiomass: ', pfba_solution.fluxes[biomass.id], '\n')
        for i in model.reactions:
            if pfba_solution.fluxes[i.id] > 1e-10 or pfba_solution.fluxes[i.id] < -1e-10:
                non_zero.append(i)
        return non_zero
    except Infeasible as i:
        print("""pFBA is infeasible, control if the coefficients of the 
        reaction equation are correct (or use a different boudary 
        reaction of the target as model objective): """, i)
         
    

def reaction_list_pruning(input_file, model, r_list):
    """
    Pruning the list of reaction with non-zero flux

    Gets the minimal set of reactions by pruning the 
    list obtained from get_raw_reaction_list function. It does so
    by excluding from the list all the exchange, transporter, cofactor
    balancing, boundary, and ATP generating reactions.

    If the pathway is suppoesd to happen in aerobiosis, the anaerobic
    specific reactions are also eliminated.
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    r_list--list of cobra.Reaction objects returned by 
        get_raw_reaction_list funciton
    Return:
    r_list:list of cobra.Reaction objects, pruned. Or Nonetype object
    """
    aerobic_to_exclude_id = ['redp', 'THD2pp', 'CYTBO3_4pp', 
                            'NADH16pp', 'EX_', 'endless', 
                            'tex', 'tpp', 'ATPS4rpp', 'DM_', 
                            'SK_', 'NTPTP', 'NADHHS', 'NADHXE']
    aerobic_to_exclude_name = ['ranspor', 'equilibr', 'anaer', 'maintenance' ]
    
    anaer_to_exclude_id = ['redp', 'EX_', 'NADK', 'FLDR2', 'endless', 'tex', 
                            'tpp', 'ATPS4rpp', 'DM_', 'SK_', 'NTPTP', 
                            'NADHHS', 'NADHXE']
    anaer_to_exclude_name = ['ranspor', 'equilibr', 'maintenance' ]
    metabolism_type = defining_metabolism(input_file, model)
    if r_list != None and r_list !=0:
        print('Raw list length = ', len(r_list))
        if metabolism_type == 'aerobic':
            for ele in aerobic_to_exclude_id:
                for i in list(r_list):
                    if ele in i.id:
                        #print(i.id, i.name, i.reaction, pfba_solution.fluxes[i.id])
                        r_list.remove(i)
            #print(len(r_list))
            
            for ele in aerobic_to_exclude_name:
                for i in list(r_list):
                    if ele in i.name:
                        #print(i.id, i.name, i.reaction, pfba_solution.fluxes[i.id])
                        r_list.remove(i)
            print('Pruned list length = ', len(r_list))
            
        elif metabolism_type == 'anaerobic':
            for ele in anaer_to_exclude_id:
                for i in list(r_list):
                    if ele in i.id:
                        #print(i.id, i.name, i.reaction, pfba_solution.fluxes[i.id])
                        r_list.remove(i)
            #print(len(r_list))
            
            for ele in anaer_to_exclude_name:
                for i in list(r_list):
                    if ele in i.name:
                        #print(i.id, i.name, i.reaction, pfba_solution.fluxes[i.id])
                        r_list.remove(i)
            print('Pruned list length = ', len(r_list))
    else:
        r_list = None
        
    return r_list


def whole_procedure_path_definition(input_file, model, raw_tsv_file_name, input_mdf_filename):
    """
    Whole pFBA analysis to get the minimal set of production reactions

    Calls the functions involved in finding the minimal set of 
    reactions for the production of the target.
    1. adds the free balancing reactions
    2. remove maintenance
    3. constrains substrate and product using the stoichiometric
        coefficients
    4. checks the carbo source to be the expected one
    5. set the target as model objective
    6. checks the metabolism
    7. gets reactions list from pFBA 
    8. processes the raw list to get the minimal reaction set
    9. generates SBtab file with the results ready for MDF calculation
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    raw_tsv_file_name--str of the intermediate .tsv SBtab file
    input_mdf_filename--str of the .tsv SBtab file for mdf analysis
    Return:
    new_r_list:list of cobra.Reaction objects, pruned, returned by 
        reaction_list_pruning funciton
    """
    substrate = get_ex_c_source_metab(input_file, model)
    target = get_ex_target(input_file, model)

    old_bounds = []
    old_bounds.append(substrate.bounds)
    old_bounds.append(target.bounds)

    # Add free balancing and remove incompatible reactions
    add_free_balancing_reactions(input_file, model)
    
    # Remove mainenance
    for i in model.reactions:
        if 'maintenance' in i.name:
            #print(i.id, i.name)
            model.remove_reactions([i.id])
            
    # Set constraints substrate
    set_constraints_c_source(input_file, model)
    print("Substrate's exchange reaction bounds : ", substrate.bounds)
    # Set constraints porduct
    set_constraints_targets(input_file, model)
    print("Target's exchange reaction bounds : ", target.bounds)
    #Check that the c-source is the expected one
    active_c_sources = cobra.flux_analysis.phenotype_phase_plane.find_carbon_sources(model)

    for index, r in enumerate(active_c_sources):
        if len(active_c_sources)==1 and substrate.id == r.id:
            continue
        elif len(active_c_sources)>1 and substrate.id == r.id:
            continue
        elif len(active_c_sources)>1 and substrate.id != r.id:
            r.lower_bound = 0
    
    print('Carbon source: ', cobra.flux_analysis.phenotype_phase_plane.find_carbon_sources(model)) 
    
    #setting the objective
    model.objective = target
    #print(model.objective.expression)
    #print(cobra.flux_analysis.phenotype_phase_plane.find_carbon_sources(model))
      
    # Eliminate the oxygen sensitive reactions
    metab_type = defining_metabolism(input_file, model)
    if metab_type == 'aerobic':            
        oxygen_sensitivity_check(model)
    elif metab_type == 'anaerobic':
        o2_conditions = model.reactions.get_by_id('EX_o2_e')
        o2_conditions.lower_bound = 0
    #por = model.reactions.POR5
    #trd = model.reactions.TRDR
    #print(por.bounds, trd.bounds)
    #print(cobra.flux_analysis.phenotype_phase_plane.find_carbon_sources(model))    
    
    # Get reactions list from pFBA and prune it
    r_list = get_raw_reaction_list(input_file, model)
    new_r_list = reaction_list_pruning(input_file, model, r_list)

    substrate.bounds = old_bounds[0]
    target.bounds = old_bounds[1]

    biomass = get_biomass_equation(model)
    model.objective = biomass

    if new_r_list != None:
        generate_SB_tab(input_file, model, new_r_list, raw_tsv_file_name, input_mdf_filename)
    else:
        print('The thermodynamic analysis cannot proceed because of infeasible pFBA')

    return new_r_list


def compound_list(pruned_r_list):
    """
    Gets the metabolites involved in the miminal set of reactions

    ------------------------------------------------------------------
    Arguments:
    pruned_r_list:list of cobra.Reaction objects, pruned, returned by 
        reaction_list_pruning funciton
    
    Return:
    compounds1: list of cobra.Metabolite objects
    """
    compounds1 = []
    for i in pruned_r_list:
        for k in i.metabolites:
            if k not in compounds1:
                compounds1.append(k)
    return compounds1


def get_c_source_concentrations(input_file):
    """
    Set range concentration range for the carbon source

    The concentration range of the metabolites are 
    needed for MDF calculations. Those have to be added in 
    the input file under c_source_concentration_range
    ------------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like

    Return:
    c_range: dictionary with 'min' and 'max' as keys and the 
        intracellular concentrations in uM as values
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            c_range = {}
            for row in reader:
                if row['c_source_concentration_range'] != '':
                    concentrations = str(row['c_source_concentration_range']).split()
                    c_range['min'] = concentrations[0]
                    c_range['max'] = concentrations[1]
    return c_range


def write_tsv_file(input_file, raw_tsv_file_name, pruned_r_list, compounds, essential_dict):
    """
    write .tsv file with the intermediate SBtab with the minimal raction set

    Write the information on reactions, metabolites, concentrations in
    order to creat a SBtab file combatible to the equilibrator-pathway
    package https://gitlab.com/equilibrator/equilibrator-pathway/-/tree/develop
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    raw_tsv_file_name--str of the intermediate .tsv SBtab file
    pruned_r_list:list of cobra.Reaction objects, pruned, returned by 
        reaction_list_pruning funciton
    compounds--list of cobra.Metabolite objects returned by the function
        compunds_list
    essential_dict--dict derived from pruned_r_list
    
    Return:
    .tsv file with the name corresponding to raw_tsv_file_name string
    """
    out_path = raw_tsv_file_name
    with open(out_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            #SBtab Reaction
            writer.writerow(["!!SBtab", "TableID='Reaction'", "TableType='Reaction'"])
            writer.writerow(["!ID", "!ReactionFormula"])
            for i in pruned_r_list:
                react = i.reaction

                if ' --> ' in react:
                    new_react = react.replace(' --> ', ' <=> ')
                    writer.writerow([i.id.lower(), new_react])
                elif ' <-- ' in react:
                    new_react = react.replace(' <-- ', ' <=> ')
                    writer.writerow([i.id.lower(), new_react])
                elif ' <=> ' in react:
                    writer.writerow([i.id.lower(), react])

            writer.writerow([])

            #SBtab Compounds
            writer.writerow(["!!SBtab", "TableID='Compound'", "TableType='Compound'"])
            writer.writerow(["!ID", "!Identifiers"])
            for i in compounds:

                if 'kegg.compound' in i.annotation.keys():
                    identifier = i.annotation['kegg.compound']
                    if type(identifier)!=list:
                        writer.writerow([i.id, "kegg:{}".format(identifier)])
                    else:
                        writer.writerow([i.id, "kegg:{}".format(identifier[0])])
                else:
                     "MissingKEGG_ID: Compound {} does not have any assigned KEGG Identifier".format(i)

            writer.writerow([])

            #SBtab Flux
            writer.writerow(["!!SBtab", "TableID='Flux'", "TableType='Quantity'", ":Unit='mM/s'"])
            writer.writerow(["!QuantityType", "!Reaction", "!Value"])
            for key, val in essential_dict.items():
                writer.writerow(['rate of reaction', key.lower(), val])

            writer.writerow([])

            #SBtab ConcentrationConstraint
            writer.writerow(["!!SBtab", "TableID='ConcentrationConstraint'", "TableType='Quantity'", ":Unit='mM'"])
            writer.writerow(["!QuantityType", "!Compound", "!Min", "!Max"])

            conc_dict = get_c_source_concentrations(input_file)
            c_source = get_c_source(input_file)

            for i in compounds:
                if c_source in i.id:
                    writer.writerow(['concentration', i.id, float(conc_dict['min']), float(conc_dict['max'])])
                elif i.id == 'pi_c':
                    writer.writerow(['concentration', i.id, 10, 10])
                elif i.id == 'ppi_c':
                    writer.writerow(['concentration', i.id, 1, 1])
                elif i.id == 'coa_c':
                    writer.writerow(['concentration', i.id, 0.1, 5])
                elif i.id == 'atp_c':
                    writer.writerow(['concentration', i.id, 5, 5])
                elif i.id == 'adp_c':
                    writer.writerow(['concentration', i.id, 0.5, 0.5])
                elif i.id == 'amp_c':
                    writer.writerow(['concentration', i.id, 0.5, 0.5])
                else:
                    writer.writerow(['concentration', i.id, 0.001, 10])


def config_SB_tab(raw_tsv_file_name, input_mdf_filename):
    """
    Adds configuration info needed for the SBdocument for MDF analysis

    It writes the fixed configuration information to the SBtab file 
    generated by write_tsv_file function.
    ------------------------------------------------------------------
    Arguments:
    raw_tsv_file_name--str of the intermediate .tsv SBtab file
    input_mdf_filename--str of the .tsv SBtab file for mdf analysis
    """
    Sd_config = SBtab.read_csv('configuration.tsv', 'Sd_config', xlsx=False)
    tab_config = Sd_config.get_sbtab_by_id('Configuration')

    # create an SBtab Document Object Sd
    Sd = SBtab.SBtabDocument('for_MDF')
    Sd.set_name('for_MDF.tsv')
    Sd.add_sbtab(tab_config)

    # you can add further SBtab tables or strings to the document
    # by using the functions add_sbtab() or add_sbtab_string()
    
    out_path = raw_tsv_file_name

    Sd_full1 = SBtab.read_csv(out_path, 'path', xlsx=False)

    tab_rea = Sd_full1.get_sbtab_by_id('Reaction')
    Sd.add_sbtab(tab_rea)

    tab_com = Sd_full1.get_sbtab_by_id('Compound')
    Sd.add_sbtab(tab_com)

    tab_flux = Sd_full1.get_sbtab_by_id('Flux')
    Sd.add_sbtab(tab_flux)

    tab_cc = Sd_full1.get_sbtab_by_id('ConcentrationConstraint')
    Sd.add_sbtab(tab_cc)

    mdf_path = input_mdf_filename
    Sd.write(filename = mdf_path)


def generate_SB_tab(input_file, model, pruned_r_list, raw_tsv_file_name, input_mdf_filename):
    """
    Generate the input file for MDF calculation

    It uses the minimal reaction set to write the tsv file used as
    input for MDF calculation. 
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    pruned_r_list:list of cobra.Reaction objects, pruned, returned by 
        reaction_list_pruning funciton
    raw_tsv_file_name--str of the intermediate .tsv SBtab file
    input_mdf_filename--str of the .tsv SBtab file for mdf analysis
    """
    compounds = compound_list(pruned_r_list) 
    #for metab in compounds:
        #if 'kegg.compound' not in metab.annotation.keys():
            #print('\nThis coupound {} does not exist in equilibrator cache'.format(metab.id))
            #print("""\nTo solve this issue you may want to look in kegg website for the id 
            # correspondant to {}, then add it as values to 'kegg.compound' 
            # key right after loading the model""".format(metab.name))
            #TODO This long message could become a warning
        #elif 'inchi_key' not in metab.annotation.keys():
            #print('This coupound {} does not have inchi_key'.format(metab.id))
    
    for i in list(pruned_r_list):
        out = i.check_mass_balance()
        if out != {}:
            print('\nThe reactions {} is not balanced'.format(i.id))
    
    fba = model.optimize()
    essential_dict={}
    for i in pruned_r_list:
        essential_dict[i.id] = fba.fluxes[i.id]
        
    write_tsv_file(input_file, raw_tsv_file_name, pruned_r_list, compounds, essential_dict)
    config_SB_tab(raw_tsv_file_name, input_mdf_filename)


def get_mdf_value(input_mdf_file):
    """
    Calculate MDF

    Uses equilibrator-pathways function to calculate MDF on 
    the input SBtab with the information on the pathway. It tries
    to run the calculation and if exceptions arise the value of 
    MDF is set to None.

    Exceptions:
    Infeasible: pFBA could be infeasible
    ParseException: there might be compounds without KEGG IDs
    MissingDissociationConstantsException: some coumpounds might not
       be associated to InChI keys. 
    ------------------------------------------------------------------
    Arguments:
    input_mdf_file--str of the .tsv SBtab file for mdf analysis

    Return:
    mdf_value: mdf value (float) and its associated magnitude or None
    """
    mdf_path = input_mdf_file

    comp_contrib = ComponentContribution()
    try:
        pp = Pathway.from_sbtab(mdf_path, comp_contrib=comp_contrib)
        pp.update_standard_dgs()
        mdf_sol = pp.calc_mdf()
        mdf_value = mdf_sol.mdf
    except Infeasible as i:
        print('''pFBA is infeasible, control if the coefficients of the 
        reaction equation are correct : ''', i)
        mdf_value = None    
    except ParseException as ex:
        print('There might be compunds not associated to any KEGG identifier: ',ex)
        mdf_value = None
    except MissingDissociationConstantsException as ex:
        print('There might be a compund without InChI key: ', ex)
        mdf_value = None
        
    return mdf_value

def mdf_analysis(input_file, model, raw_tsv_file_name, input_mdf_filename):
    """
    Combines all the functions for identifying the pathway and calculating MDF

    1.  calls the function for the identification fo the minimal
    set of reactions active in the conversion. 
    2. uses the gneerated SBtab for mdf calculation

    If anthing in the process goes wrong and mdf can't be calculated
    or the pathway is not identified the values for pathway length and
    mdf are set to None. 
    ------------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    raw_tsv_file_name--str of the intermediate .tsv SBtab file
    input_mdf_filename--str of the .tsv SBtab file for mdf analysis
    
    Returns:
    pathways_MDF:mdf value (float) and its associated magnitude returned
        by get_mdf_value or NoneType
    path_len: int, length of the list with the minimal reaction set as
        obtained from whole_procedure_path_definition function
    """
    pruned_r_list = whole_procedure_path_definition(input_file, model, raw_tsv_file_name, input_mdf_filename)
    if pruned_r_list != None and pruned_r_list != 0:
        pathways_MDF = get_mdf_value(input_mdf_filename)
        path_len = len(pruned_r_list)
    else:
        print('The thermodynimac analysis has been unsuccesful')
        pathways_MDF = None
        path_len = None
    return pathways_MDF, path_len
    
