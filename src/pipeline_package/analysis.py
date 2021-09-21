#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module to perform the Metabolic analysis for reaction knock-ins.

This module has the functions performing GapFilling analysis to find
knock-in reactions for allowing growth on the selected substrates
and production of target compounds as indicated in the input file. 
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"

import logging
import sys
from copy import deepcopy

import cobra
from cobra.flux_analysis import phenotype_phase_plane
from cobra.flux_analysis import gapfill
import csv

from pipeline_package.input_parser import parser, get_metabolites, get_reactions, set_bounds_ex
from pipeline_package.path_definition_mdf import * 
from pipeline_package.import_models import get_ID_reference_model, get_expression_host


logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
#import warnings

#warnings.filterwarnings('error', category=UserWarning, module='analysis.py')

growth_on_preferred_c = {}

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

def shut_down_c_sources(input_file, model):
    """
    Set uptake of default carbon sources to 0 

    Given a model 
    1) It defines what the carbon sources are with the function 
    phenotype_phase_plane.find_carbon_sources from cobra.
    2) Then it shuts down their consumption setting the lb to 0
    3) Then sets the lower bound of the exchange reaction
    of the selected carbon source to -1000
    --------------------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    """
    c_sources = phenotype_phase_plane.find_carbon_sources(model)
    print('\n---Carbon sources---')
    for n in range(len(c_sources)):
        carbon = c_sources[n]
        print(carbon.name) 
        print(carbon.id)
        print(carbon.reaction)
        print("\nOld bounds: ", carbon.bounds)
        carbon.lower_bound = 0
        print("New bounds: ", carbon.bounds, '\n')
    ex = get_ex_c_source_metab(input_file, model)
    ex.lower_bound = -1000


def iter_number(input_file):
    """
    Extracts the user-decision on the number of iteration 
    
    The user must indicate the number of iterations of GapFilling analysis
    by filling the number in correspondance of "iterations" column 
    --------------------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like

    Return:
    n: int number of iterations
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            for row in reader:
                if row['iterations'] != '':
                    n = int(row['iterations'])
                    #print('Iterations: ',n)
                    break
    return n

def get_max_growth_on_preferred_c(input_file, model):
    """
    Define the maximal growth on the preferred carbon source

    Read the input file to get the information on the preferrd carbon 
    source: the metabolite and the lower bound of the uptake reaction
    They must both be indicated by the user in the input file in 
    correspondance preferred_c_source and max_average_uptake columns.
    FBA is performed to get the flux through biomass reaction.
    --------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace

    Return:
    ub_bomass: float flux value through biomass reaction
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        biomass = get_biomass_equation(model) 
        for row in reader:
            if row['preferred_c_source'] != ''and row['max_average_uptake'] != '':
                pref_c_source = row['preferred_c_source'] 
                lb_preferred_uptake = float((row['max_average_uptake'])) 
        uptake_reaction = model.reactions.get_by_id('EX_'+pref_c_source+'_e')
        logging.debug(uptake_reaction)
        uptake_reaction.lower_bound = lb_preferred_uptake
        logging.debug(lb_preferred_uptake)
        model.solver = 'glpk'
        model.objective = biomass
        logging.debug(type(model.solver))
        fba_pref_source = model.optimize()
        for i in model.reactions:
            if i.flux <= -0.5 and "EX_" in i.id:
                print(i.id, i.reaction, i.flux, '\n')
        ub_biomass = fba_pref_source.objective_value
        logging.debug(ub_biomass)
        uptake_reaction.lower_bound = 0
        growth_on_preferred_c['G_on_preferred_c'] = ub_biomass
        logging.debug(growth_on_preferred_c)
    return ub_biomass


def set_bounds_obj(input_file, model, universal):
    """
    Set upper bound of biomass reaction and the lb of the uptake reaction for optimization of consumption

    It uses the float value resulting from get_max_growth_on_preferred_c
    function as upper bounds of biomass. Growth is then constrained to 
    a maximum that equals the upper bound, hence the consumption rate
    of the desired compound can be recorded. Additionally, the lb of the 
    exchange reaction of the selected carbon source is set to -1000.
    --------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace
    """
    biomass = get_biomass_equation(model) 
    metabolites = get_metabolites(input_file)
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        #metabolites = get_metabolites(input_file)
        
        print('Old biomass (objective) bounds = ', biomass.bounds)
        biomass.upper_bound = get_max_growth_on_preferred_c(input_file, model) # originally 0.87 in Ecoli, because it is unlikely to grow more than E.coli does on a middle-high glucose consentration
        
        #The following are not useful, since those setting are valid only for the optimization of consumption and in this case the only constraint needed is the one on the ub of biomass.

        # for metab in metabolites:
        #    for row in reader:
        #        if metab in row['metabolite_BiGG_ID'] and row['carbon_source']=='yes' and row['consumption']=='yes':
        #            biomass.upper_bound = get_max_growth_on_preferred_c(input_file, model) # originally 0.87 in Ecoli, because it is unlikely to grow more than E.coli does on a middle-high glucose consentration
        #        elif metab in row['metabolite_BiGG_ID'] and row['carbon_source']=='yes' and row['consumption']!='yes' and row['consumption']!='no':
        #            continue 
        print('New biomass (objective) bounds = ', biomass.bounds)
        
        for i in metabolites:

            if ('EX_' + i + '_e') in model.reactions:
                print('\n{} is in the medium'.format(i), '\n')
            elif ('EX_' + i + '_e') not in model.reactions:
                print('\n{} is not in the medium'.format(i), '\n')
                if ('EX_' + i + '_e') in universal.reactions:
                    exchange = universal.reactions.get_by_id('EX_' + i + '_e')
                    exchange.lower_bound = -1000
                    exchange.upper_bound = 0
                    model.add_reaction(exchange)
                else:
                    extracellular = Metabolite(i + '_e')
                    extracellular = Metabolite(i + '_e', name='extracellular' + i, compartment='e')
                    exch = Reaction('EX_' + i + '_e')
                    exch.name = 'Exchange {}'.format(extracellular)
                    exch.lower_bound = -1000
                    exch.upper_bound = 0
                    exch.add_metabolites({extracellular: -1})
                    model.add_reaction(exch)
                    print('Ther reaction {} has been added to the reference model, hence {} is now in the medium'.format(
                            'EX_' + i + '_e', i + '_e'))
            set_bounds_ex(input_file, model, i)

def iter_gf(input_file, model, universal):
    """
    Run gapfilling on model using universal as reaction database.

    This function runs gapfilling on the model using the universal of
    choice as database of reactions and it performes the number of 
    iteration defined in the input file.

    Notes: If the optimisation of the objective model is infeasible 
    the flux of the objective equals 0.

    If the optimization of the model with the user's settings is
    feasible it creates an analogous dictionary in which  
    solutions instrad of a dictionary is a float corresponding to
    the flux though biomass reaction.
    ------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace
            
    Return:
    solutions: dictionary type, where the number of the keys equals the 
        number of iterations and the values correspond to the identified
        reactions to fill the gaps, or float type

    """
    n = iter_number(input_file)
    #try:
        #fba_model_obj = model.optimize()
    #except UserWarning: 
    model.solver = 'glpk'
    logging.debug(model.objective.expression)
    #glc = model.reactions.EX_glc__D_e
    #logging.debug(glc.bounds)
    fba_consumption = model.optimize()
    logging.debug('passed here too')
    logging.debug(fba_consumption.objective_value)
    if fba_consumption.objective_value == None or fba_consumption.objective_value <= 0:    
        model.solver = 'glpk'
        logging.debug(type(model.solver))
        result = gapfill(model, universal, demand_reactions=False, exchange_reactions=False, iterations=n) 
        solutions ={}
        for i, entries in enumerate(result):
            reacts = []
            for e in entries:
                reacts.append(e.id)
                if "glc__D_c" in e.reaction:
                    print("Reaction ID: ", e.id, "But it might be a circulation")
                    solutions["Run {}".format(str(i + 1))] = reacts
                else:
                    solutions["Run {}".format(str(i + 1))] = reacts
    else:
        print('The model can already satisfy the objective')
        solutions = fba_consumption.objective_value
    return solutions
    #TODO: handle errors, what if no solution found? e.g None
    
def get_c_source(input_file):
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel')
            c_source = ''
            for row in reader:
                if row['carbon_source'] == 'yes':
                    c_source += row['metabolite_BiGG_ID']
    return c_source

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

def get_ID_non_c_source_metabs(input_file, model):
    """
    List the exchange reaction of the metabolites that are not c-sources

    This function gets the IDs of the metabolites that have to be consumed
    but are not carbon sources. It does not account for the target.
    --------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace

    Return:
    non-carb: list of ID of non-c-source metabolites to be consumed.
    """
    to_consume = consumption_metabs(input_file)
    logging.debug(to_consume)
    non_carb =[]
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        for row in reader:
            for i in to_consume:
                if i in row['metabolite_BiGG_ID'] and row['carbon_source'] == 'no':
                    non_carb.append(i)
                else:
                    pass
    return non_carb

def consumption_metabs(input_file):
    """
    List of metabolites that are consumed.

    It includes carbon containing metabolites as well as non c-source
    ones, as long as in the input file it is indicated that they 
    should be consumed.
    --------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like

    Return:
    consumption: list of IDs of metabolite to be consumed.
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
                reader = csv.DictReader(csvfile, dialect='excel')
                metabs = get_metabolites(input_file)
                consumption = []
                for row in reader:
                    for i in metabs:
                        if (i in row['metabolite_BiGG_ID'] and 
                                row['consumption'] in ('yes', 'Yes', 'Y', 'y', 'YES')):
                            if i not in consumption:
                                consumption.append(i)
                        elif (i in row['metabolite_BiGG_ID'] and  # for accounting numbers in cell
                                row['consumption'] != 'no' and 
                                row['consumption'] != ''):
                            if i not in consumption:
                                consumption.append(i)

    return consumption

def remove_rlist(l_reacts, model):
    """
    Eliminate a list of reactions from the model.

    --------------------------------------------
    Arguments:
    l_reacts--list containing the BiGG IDs of the reactions to remove
    model--cobra.Model reference model in BiGG namespace
    """
    for i in l_reacts:
        model.remove_reactions([i])

def eval_sol(input_file, model, universal):
    """
    Optimize for consumption and evaluate the solutions

    1) Performs GF
    2) evaluates the output by creating a dictionary with
        - number of iteration
        - models
        - reactions added
        - FBA for growth
        - flux through added reactions
        - flux through exchange reaction of metabolites to consume
        - flux through exchange reaction of metabolites to produce
    --------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace

    Return:
    output: dict with the information on each solution for 
        optimized consumption
    
    """
    biomass = get_biomass_equation(model)
    model.objective = biomass
    logging.debug(model.objective.expression)
    logging.debug(biomass.bounds)
    #right_model = deepcopy(model)
    iterations = iter_number(input_file) 
    logging.debug('passed here')
    print('Starting reaction search with GapFilling . . .')
    solutions = iter_gf(input_file, model, universal)
    metabs = get_metabolites(input_file)
    reacts = get_reactions(input_file)
    #substitute the following two lines with consumption_metabs and get_production_objectives
    #ex_c_source = get_ex_c_source_metab(input_file, model)
    #ex_non_c_metabs = get_ex_non_c_source_metabs(input_file, model)
    to_consume = consumption_metabs(input_file)
    to_produce = get_production_objectives(input_file)
    output = {}
    models={}
    for n in range(1,iterations+1): #number of iterations
        consume ={}
        produce ={}
        models={}
        added_reacts={}
        print('\n---Model {}---'.format(n))
        logging.debug(type(model.solver))
        models['model']='Model{}'.format(n)
        if type(solutions)==dict and len(solutions.keys())!=0: #runs the analysis of the results only if reaction finding was needed.
            sol = solutions['Run {}'.format(n)]
            m = len(sol) #number of reactions as solution fo each iteration
            for i in range(m):
                    identifier = sol[i]
                    react = universal.reactions.get_by_id(identifier)
                    model.add_reaction(react)
                    print('\nReaction {}, solution of round {} has been added to the model'.format(react.id, n))
            model.solver = 'glpk'
            logging.debug(type(model.solver))
            fba_model_x = model.optimize()
            print('\nGrowth rate: ', biomass.flux)
            for x in metabs:
                r=model.reactions.get_by_id('EX_'+x+'_e')
                print('\nThe flux throughr {} is: '.format(r.id), r.flux)
            for i in range(m):
                identifier = sol[i]
                r = model.reactions.get_by_id(identifier)
                print('\nThe flux throughr {} is: '.format(identifier), r.flux)
                added_reacts[sol[i]]=r.flux
            for compound in to_consume:
                exchange = model.reactions.get_by_id('EX_'+compound+'_e')
                consume[exchange.id] = exchange.flux
            for target in to_produce:
                exch = model.reactions.get_by_id('EX_'+target+'_e')
                produce[exch.id] = exch.flux

            #carbon_flux[ex_c_source.id] = ex_c_source.flux
            #for x in ex_non_c_metabs:
            #    flux_non_c_metabs[x.id] = x.flux
            info_models = (sol, biomass.flux, added_reacts, consume, produce)
            #(reactions, growth_rate, added_reactions, carbon, non_carbon)=info_models
            output[n] = (models, info_models)
            remove_rlist(sol, model)
        else: #if the solution of iter_gf was the flux value though biomass 
                exch_c_source = get_ex_c_source_metab(input_file, model)
                fba_model = model.optimize()
                for target in to_produce:
                    exch = model.reactions.get_by_id('EX_'+target+'_e')
                output[1] = ({'model':'Model1'}, (None, solutions, None, {exch_c_source.id: exch_c_source.flux}, {exch.id: exch.flux}))
    return output

def remove_duplicated_sol(output_dict):
    """
    Remove repeated solutions for optimised consumption.

    Remove repeated solutions resulting from GapFilling analysis 
    for optimised consumption. GapFilling algorythm implemented in
    Cobra https://cobrapy.readthedocs.io/en/latest/gapfilling.html
    assigns penalyies to solutions already used but it can 
    happen that using the same solution is preferred to finding new 
    ones, e.g. involving a higher number of reactions. Hence this 
    function records only the unique solutions
    --------------------------------------------
    Argument:
    output_dict--dict returned by eval_sol function
    
    Return:
    output_dict: dict, modified solution dictionary
    """
    solutions_gf = []
    for i in list(output_dict.values()):
        if type(i[1][0]) != float and i[1][0]!=None: #output processing is needed only if there are reaction addition as solution
            if i[1][0] not in solutions_gf:
                solutions_gf.append(i[1][0])
            elif i[1][0] in solutions_gf:
                del output_dict[int(i[0]['model'][-1:])]
    return output_dict

def generate_output_file(processed_output):
    """
    Create intermediate csv file with the solutions for optimised consumption
    
    Create file with .csv extension by reading the processed dictionary 
    resulting from remove_duplicated_sol function. It contains 
    information such us model's version, added reactions and FBA results.
    --------------------------------------------
    Arguments:
    processed_output--dict returned by remove_duplicated_sol function
    out_file--str name of the intermediate .csv file 
    
    Return:
    saves a csv file with the idicated name
    """
    with open('pipeline/outputs/consumption.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'model_version', 'growth_rate', 
            'reactions_added', 'fluxes', 
            'compounds_to_consume', 'uptake_rate', 
            'compounds_to_produce', 'production_rate'
            ])
        for i in processed_output.values():
            if type(i[1][0]) != float and i[1][0]!=None:#the file is created only if there are reaction addition as solution
                cs = i[1][3]
                exc = []
                val = []
                mod = [i[0]['model']]
                for x, y in cs.items():
                    exc.append(x)
                    val.append(y)

                tar = i[1][4]
                exp = []
                valp = []
                for prod, flu in tar.items():
                    exp.append(prod)
                    valp.append(flu)
                
                r = []
                f = []
                for react, flux in i[1][2].items():
                    r.append(react)
                    f.append(flux)
                
                for n in range(len(r)):
                    if (n+1) <= len(exc) and (n+1) <= len(mod):
                        writer.writerow([(i[0]['model']), i[1][1], r[n], f[n], exc[n], val[n],exp[n], valp[n]])
                    elif (n+1) <= len(exc) and (n+1) > len(mod):
                        writer.writerow(['', '', r[n], f[n], exc[n], val[n],'',''])
                    elif (n+1) > len(exc) and (n+1) > len(mod):
                            writer.writerow(['', '', r[n], f[n], '','','',''])

def analysis_gf_sol(input_file, model, universal):
    """
    Call functions for optimizing a model for growth on unusual c-source

    The following steps are performed when calling this function:
    1) Sets the right bounds for biomass
    2) automatically calls functions for GapFilling and evaluating the results 
    3) creates the output file in csv format
    --------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace

    Return:
    sol_for_csv: dict with the information on unique solutions 
        for optimized consumption. 
    """
    metabolites = get_metabolites(input_file)
    #for metab in metabolites:
      #  r = model.reactions.get_by_id('EX_'+metab+'_e')
        #print(r.id, " | ", r.reaction, " | ", r.bounds)
    
    set_bounds_obj(input_file, model, universal) # Onyl biomass is constrained by setting ub = growth on preferred substrate
    #shut_down_c_sources(input_file, model)
    output = eval_sol(input_file, model, universal)
    sol_for_csv = remove_duplicated_sol(output)
    generate_output_file(sol_for_csv)
    return sol_for_csv

def get_production_objectives(input_file):
    """
    List the target metabolites to produce

    It reads the input file and gets which metabolites has to
    be produced.
    --------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like
    
    Return:
    production_objects: list of BiGG IDs of the target metabolites  
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        production_objects = []
        for row in reader:
            if (row['production'] in ('yes', 'Yes', 'Y', 'y', 'YES') and
                    row['metabolite_BiGG_ID']!= ''):
                production_objects.append(row['metabolite_BiGG_ID'])
    return production_objects    
    
def info_fba_optimization(input_file, model, target, sol):
    """
    Print fluxes through relevan reactions if the target can be produced

    If FBA with the exchange reaction of the target as objective
    is feasible, this function prints the reaction rates of 
    relevant reactions 
    --------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    target--str target metabolite BiGG ID
    sol--list of BiGG IDs of reactions found as solutions of
       GapFilling for optimized growth on the indicated source.
    """
    print('---{}---'.format(target.upper()))
    ex_prod = model.reactions.get_by_id('EX_'+target+'_e')
    print('\nA flux through the exchange reaction for {} is {} without the need of GapFilling'.format(target.upper(), ex_prod.flux))
    biomass = get_biomass_equation(model)
    print('\nThe growth rate is: ', biomass.flux)
    metabs = get_metabolites(input_file)
    for m in metabs:
            r=model.reactions.get_by_id('EX_'+m+'_e')
            print('\nThe flux throughr {} is: '.format(r.id), r.flux)
    for x in sol:
        r = model.reactions.get_by_id(x)
        print('\nThe flux throughr {} is: '.format(x), r.flux)

def add_sol_cons(dict_sol_consumption, model_number, model, universal):
    """
    TODO: eliminate or modify, not working
    Add reaction to allow growth on selected compound.

    This function adds to the model the reactions which are
    solutions for the optimization of consumption.
    --------------------------------------------
    Arguments:
    dict_sol_consumption--dict as outptu from analysis_gf_sol function
    model_number--int identfying iteration of GapFilling
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace
    """
    for reaction in dict_sol_consumption[model_number][1][0]:
        to_add = universal.reactions.get_by_id(reaction)
        model.add_reaction(to_add)

def production_analysis(input_file, model, universal):
    """
    Calls function for GapFilling analysis with target production as objective
    
    This function 
    1) changes the objective of the model to the exchange of the
    target compound(s)
    2) Calls Gapfilling function and creates a dict with the solutions 
    resulting from the iterations for each target compound
    --------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace

    Return:
    dict_products: dict with the reactions resulting from GapFilling
    """
    to_produce = get_production_objectives(input_file)
    dict_products={}
    for compound in to_produce:
            prod = model.reactions.get_by_id('EX_'+compound+'_e')
            model.objective = prod
            sol_for_p = iter_gf(input_file, model, universal)
            dict_products[compound]=sol_for_p
    if type(dict_products[compound]) == dict:
        return dict_products
    else: #TODO: fix, I don't get why I did it like this. Check iter_gf todo
        return dict_products[compound]

def set_constraint_production(input_file, model):
    """
    Set constraints on uptake and min growth rate
    
    Given the input file and the model it constraints uptake 
    and min growth rate:
    1) uptake:
    - if the user has chosen a value for the lb of the carbon source
        it sets the selected one 
    - if the user hasn't set any value, a message is given and 
        the lb is set default to -100 
    2) growth rate:
    - the lb of growth rate is set to be 5% of the growth on 
        the preferred substrate
    --------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace

    Return:
    model: cobra.model with the modified/constrained reactions
    """
    biomass = get_biomass_equation(model) 
    
    model.solver = 'glpk'
    to_consume = consumption_metabs(input_file)
    exch_c_source = get_ex_c_source_metab(input_file, model)
    #logging.debug(exch_c_source)
    c_exch_bounds = exch_c_source.bounds
    c_source = get_c_source(input_file)
    non_c_met=get_ID_non_c_source_metabs(input_file, model)
    #logging.debug(c_source)
    #logging.debug(non_c_met)
    old_biomass = biomass.bounds
    print('\nBounds of biomass during optimization of consumption = ', biomass.bounds)
    growth_on_p_subs = growth_on_preferred_c['G_on_preferred_c']
    biomass.lower_bound = (growth_on_p_subs/100)*5
    print('\nBounds of biomass during optimization of production = ', biomass.bounds)
    new_biomass = biomass.bounds
    
    non_c=[]
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        logging.debug('passed here')
        reader = csv.DictReader(csvfile, dialect='excel')
        
        #Set constraints
        for row in reader:
            for metab in to_consume:
                logging.debug(metab)
                if metab == c_source:
                    logging.debug(metab)
                    if metab in row['metabolite_BiGG_ID'] and row['carbon_source']=='yes' and row['consumption']=='yes':
                        lb_uptake_c = -100
                    elif metab in row['metabolite_BiGG_ID'] and row['carbon_source']=='yes' and row['consumption']!='yes' and row['consumption']!='no':
                        lb_uptake_c = float(row['consumption'])
                        logging.debug(lb_uptake_c)
                    exch_c_source.lower_bound = lb_uptake_c
                else:
                    if metab in row['metabolite_BiGG_ID'] and row['carbon_source']=='no' and row['consumption']=='yes':
                        lb_uptake = -1000
                        non_c.append(lb_uptake)
                    elif metab in row['metabolite_BiGG_ID'] and row['carbon_source']=='no' and row['consumption']!='yes' and row['consumption']!='no':
                        lb_uptake = float(row['consumption']) 
                        non_c.append(lb_uptake)
                    #print('\nBounds uptake {} = '.format(metab), exch_c_source.bounds)
            #Get the balanced equation from input file
            if row['balanced_equation'] != '':
                equation = row['balanced_equation']
        
    for ind, m in enumerate(non_c_met):
        exch = model.reactions.get_by_id('EX_'+m+'_e')
        exch.lower_bound = non_c[ind]
        
    p_exch = exch_c_source.bounds
    logging.debug(p_exch)
    #get species name of the expression host
    host_species = get_expression_host(input_file)
    #get expression host model BiGG ID
    model_ID = get_ID_reference_model(input_file)

    with open('pipeline/outputs/detailed_output.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'constraints_consumption', 'constraints_production',
            'balanced_equation', 'expression_host', 'host_model_ID' 
            ])
        writer.writerow([
            '#bounds of the biomass reactions and the exchange reaction of the substrate when doing reactions search for growth on the substrate', '#bounds of the biomass reactions and the exchange reaction of the substrate when doing reactions search for production of the target',
            '#balanced equation for the conversion of the substrate into the target product', "#name of the expression host species, ", "#BiGG identifier of the expression host's model" 
            ])
        writer.writerow([
            'Biomass_reaction={}'.format(old_biomass), 
            'Biomass_reactions={}'.format(new_biomass),
            equation, host_species, model_ID
            ])
        writer.writerow([
            'Substrate_uptake={}'.format(c_exch_bounds), 'Substrate_uptake={}'.format(p_exch)
            ])
        writer.writerow(['---------------------------------------------'])        
      

          

def dict_prod_sol(input_file, sol_for_csv, model, universal):
    """
    Call functions for optimization of production and store information 

    This function generates a dictionary with the solutions for the 
    production of the target compound. 
    For each unique solution allowing the model to grow on the 
    indicated substrate it runs GapFilling and creates a dictionary
    with the output. 
    --------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    sol_for_csv--dict as output from analysis_gf_sol function
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace

    Return:
    dict_sol_prod: dict with the solutions of GapFilling for production
        associated to each model viariant growing on the desired 
        substrete    
    """
    biomass = get_biomass_equation(model)
    set_constraint_production(input_file, model) # set max uptake as indicated by the user 
    #full_out = eval_sol(input_file, model, universal1)
    #processed_out = remove_duplicated_sol(full_out)
    to_produce = get_production_objectives(input_file)
    print('\nThe metabs to produce are: ', to_produce)
    logging.debug(type(model.solver))
    dict_sol_prod = {}
    for i in sol_for_csv:
        dict_for_models = {}
        if sol_for_csv[i][1][1] > 0:
            logging.debug(sol_for_csv[i][1][0])
            if sol_for_csv[i][1][0]!= None:
                for reaction in sol_for_csv[i][1][0]:
                    to_add = universal.reactions.get_by_id(reaction)
                    model.add_reaction(to_add)
                logging.debug('passed through here')
                for reaction in sol_for_csv[i][1][0]:
                    added = model.reactions.get_by_id(reaction)
                    logging.debug(added)
            #add_sol_cons(sol_for_csv, i, model, universal)
            model.solver = 'glpk'
            #logging.debug(type(model.solver))
            model.optimize()
            
            #if biomass.flux < 0.426:
            #    biomass.lower_bound = biomass.flux
            #elif biomass.flux > 0.426:
            #    biomass.lower_bound = 0.426
                
            print("\n---"+str(sol_for_csv[i][0]['model'][-1:])+"---")

            sol_gf_production = production_analysis(input_file, model, universal)
            
            if type(sol_gf_production) == str:
                for target in to_produce:
                    info_fba_optimization(input_file, model, target, sol_for_csv[i][1][0])
            else:
                dict_for_models[sol_for_csv[i][0]['model'][-1:]]=sol_gf_production
            if sol_for_csv[i][1][0]!= None:
                list_reactions_cons = sol_for_csv[i][1][0]
                logging.debug(list_reactions_cons)
                remove_rlist(list_reactions_cons, model)
            dict_sol_prod[sol_for_csv[i][0]['model'][-1:]] = dict_for_models
        else:
            pass
    
    return dict_sol_prod
    
    

def cons_prod_dict(input_file, model, universal, sol_cons, sol_prod):
    """
    Create a dictionary with the reactions knock-ins for both optimizations

    If the model can already grow on the selected carbon source and or 
    produce the indicate target, the dictionary will be generated with the 
    fluxes derived from the optimization without the need of the addition 
    of any reaction. 
    -----------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace
    sol_cons--dict retured by analysis_gf_sol function
    sol_prod--dict retured by dict_prod_sol function

    Return:  
    candp: dict with comprehensive information on reactions additions 
        and fluxes for each model veraint growing on the indicated
        substrate and producing the target.

        Keys = consumption_n where n = number in range(len(sol_cons))
               production_n where n = number in range(len(sol_cons))
                                These keys are both dictionaries.
    """
    biomass = get_biomass_equation(model)
    to_consume = consumption_metabs(input_file)
    to_produce = get_production_objectives(input_file)
    candp = {}
    logging.debug(type(model.solver))
    model.solver = 'glpk'
    logging.debug(type(model.solver))
    m=0
    for i in sol_cons:
        m+=1
        if sol_cons[int(i)][1][0]!=None:
            #print('------------ consumption sol', i, '--------------')
            #add reactions for consumption
            for reaction in sol_cons[int(i)][1][0]:
                    #print('\n', reaction)
                    to_add = universal.reactions.get_by_id(reaction)
                    model.add_reaction(to_add)
            #add_sol_cons(sol_cons, i, model, universal)
                    logging.debug(model.reactions.get_by_id(reaction))
            #print('\n-----production sol ', i, '-----')
            # differentiate between targets  
            #mright_reactions = deepcopy(model)      
        dict_for_target={}
        if type(sol_prod[str(i)][str(i)])!= float:
            for target in sol_prod[str(i)][str(i)].keys():
                #print('\n target = ', target)
                prod = model.reactions.get_by_id('EX_'+target+'_e')
                
                
                #clean out the redundant gf solutions for production
                gfprod = []
                for x in sol_prod[str(i)][str(i)][target].values():
                    if x not in gfprod:
                        gfprod.append(x)
                
                #add reactions for production 
                for n in range(len(gfprod)):
                    #print('\n reaction for production', gfprod[n])
                    for l in range(len(gfprod[n])):
                        #print(gfprod[n][l])
                        #for reaction in j:
                        #print('\n reactions for production added', gfprod[n][l])
                        toadd = universal.reactions.get_by_id(gfprod[n][l])
                        model.add_reaction(toadd)
                    model.objective = prod
                    fba_all = model.optimize()
                    #print('\n', model.objective.expression, fba_all.objective_value)
                    consumption={}
                    for metab in to_consume:
                        exchange = model.reactions.get_by_id('EX_'+metab+'_e')
                        consumption['EX_'+metab+'_e'] = exchange.flux
                    production={}
                    for t in to_produce:
                        ex = model.reactions.get_by_id('EX_'+t+'_e')
                        production['EX_'+t+'_e'] = ex.flux
                    added_reactions={}
                    if sol_cons[int(i)][1][0]!=None:
                        for solution in sol_cons[int(i)][1][0]:
                            print('\n solutions for consumption = ', solution)
                            solc = model.reactions.get_by_id(solution)
                            added_reactions[solution]=solc.flux
                    for j in gfprod[n]:
                        #for reaction in j:
                            solp = model.reactions.get_by_id(j)
                            added_reactions[j]=solp.flux
                    
                    # thermodynamic analysis
                    # make a copy of the model
                    mthermo = deepcopy(model) 
                    mdf_value, path_length = mdf_analysis(input_file, mthermo, './pipeline/outputs/{}_Run_{}.tsv'.format(target, str(n+1)), './pipeline/outputs/input_mdf_{}_Run_{}.tsv'.format(target, str(n+1)))
                    thermodynamic = {} 
                    if mdf_value == None:
                        thermodynamic['mdf'] = mdf_value
                    else:
                        thermodynamic['mdf'] = mdf_value.magnitude
                    thermodynamic['pathway_length'] = path_length

                    # make a copy of the model and change setting to simulation uptake
                    mcopy = deepcopy(model) 
                    # change objective to growth:
                    bio_copy = get_biomass_equation(mcopy)
                    mcopy.objective = bio_copy
                    # biomass upper bound 
                    bio_copy.upper_bound = growth_on_preferred_c['G_on_preferred_c']
                    # no constraint in uptake 
                    ex_c_source = 'EX_' + get_c_source(input_file) + '_e'
                    r_ex = mcopy.reactions.get_by_id(ex_c_source) 
                    r_ex.lower_bound = -1000
                    # FBA and record the flux through the exchange reaction of the carbon source
                    fba_consumption_2 = mcopy.optimize()
                    logging.debug('objective value mcopy')
                    logging.debug(fba_consumption_2.objective_value)
                    #object with final consumption value
                    consumption_final = {}
                    for metab in to_consume:
                        exchange = mcopy.reactions.get_by_id('EX_'+metab+'_e')
                        consumption_final['EX_'+metab+'_e'] = exchange.flux
                    #consumption_final = fba_consumption_2.fluxes[ex_c_source]

                    dict_for_target[target+'_'+'Run_'+str(n+1)]=(gfprod[n], fba_all.fluxes[biomass.id], consumption, production, added_reactions, thermodynamic, consumption_final, mcopy)
                    list_reactions_prod = gfprod[n]
                    remove_rlist(list_reactions_prod, model)
                    logging.debug('end round')
            
            candp['consumption_'+str(i)] = (sol_cons[i]) 
            candp['production_'+str(i)] = dict_for_target

        else:
            candp['consumption_'+str(m)] = (sol_cons[i])
            # thermodynamic analysis
            # make a copy of the model
            mthermo = deepcopy(model) 
            for target in to_produce:
                mdf_value, path_length = mdf_analysis(input_file, mthermo, './pipeline/outputs/{}_formdf.tsv'.format(target), './pipeline/outputs/final_input_mdf_{}.tsv'.format(target))
                thermodynamic = {} 
                if mdf_value == None:
                    thermodynamic['mdf'] = mdf_value
                else:
                    thermodynamic['mdf'] = mdf_value.magnitude
                thermodynamic['pathway_length'] = path_length 
                candp['production_'+str(m)] = {target:{'EX_'+to_produce[0]+'_e flux': sol_prod[str(i)][str(i)], 'thermodynamic':thermodynamic}}  

        if sol_cons[int(i)][1][0]!=None:
            list_cons = sol_cons[int(i)][1][0]
            remove_rlist(list_cons, model)
    return candp

   