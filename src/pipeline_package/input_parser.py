#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module to read csv input file and initialize the models.

This module reads the information from the input file and it 
initialize the reference and universal model accordingly
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"

import cobra
from cobra import Model, Reaction, Metabolite
import csv
import sys
import logging
logging.basicConfig(stream=sys.stderr, level=logging.ERROR)


def get_metabolites(input_file):
    """
    Get the metabolite contained in the input file in CSV format

    Input file must be in the same folder as the script. 
    The file name must be a string including .csv extension.
    --------------------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like

    Return:
    metabolites: list of metabolites IDs as written in the input file
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        metabolites = []
        for row in reader:
            if row['metabolite_BiGG_ID'] != '':
                x = row['metabolite_BiGG_ID']
                if '#' not in x:
                    metabolites.append(x)

    return metabolites

def get_reactions(input_file):
    """
    Get all the reactions from the input file

    Reads the input file and gets the reactions that the user might 
    want to add to the universal model. 
    Input file must be in the same folder as the script. 
    The file name must be a string including .csv extension.
    --------------------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like

    Return:
    react_IDs: list of BiGG IDs of the reactions to be added to the
        universal model
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        react_IDs = []
        for row in reader:
            if row['Reaction_ID'] != '':
                r = row['Reaction_ID']
                if '#' not in r:
                    react_IDs.append(r)
    return react_IDs

def check_if_metabs_in_universal(input_file): 
    """
    Check if the metabolites in input file are all in universal model 

    Check if the reactions in the input file contains metabolites 
    that are not included in the universal model.
    --------------------------------------------------------
    Argument:
    input_file--str input file in .csv format dictionary like

    Return:
    react_w_non_bigg_metabs: list of BiGG IDs of the reactions to be 
        added to the universal model containing metabolites absent in
        BiGG.
    """
    # TODO: delete since unused
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        react_w_non_bigg_metabs = []
        for row in reader:
                if row['metabolite_name'] == '': # Which means is empty
                    # print('All metabolites exist in the universal DB for reaction {}'.format(row['Reaction_ID']))
                    continue
                elif row['metabolite_name'] != '' and '#' not in row['metabolite_name']:
                    react_w_non_bigg_metabs.append(row['Reaction_ID'])
                    #print('The reaction {} contains metbolites that had to be added to the universal model'.format(row['Reaction_ID']))
    return react_w_non_bigg_metabs

def add_transport_with_metabs_already_in_universal(input_file, metab, universal, model):
    """
    Add transport reaction to reference model. 
    
    The reaction must be written in a whole line of the CSV file.
    The reactants and products must be separated by a space and their
    coefficients must be in the same order.
    To write numbers in a cell of excell in a non-equation format an 
    apex ' has to be added at the beginning of the cell.
    -------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    metab--str ID of the metabolite
    universal--cobra.Model universal model
    model--cobra.Model reference model
    """
    #TODO: delete since unused
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        metabolites = get_metabolites(input_file)
        if metab in metabolites:
            for row in reader:
                if ("transp" in row['Reaction_name'] and 
                    metab in row['Reaction_equation']):  # TODO: FIX/IMPROVE!
                    # e.g.
                    # if (metab+'e') in row['Reaction_equation'] and (metab+'c') in row['Reaction_equation']:
                    # or [(metab+'p') in row['Reaction_equation'] and (metab+'c') in row['Reaction_equation']]:

                    # Check if any reaction has metabolite to be added to the universal
                    in_uni = check_if_metabs_in_universal(input_file)
                    if row['Reaction_ID'] not in in_uni:

                        # Reading the metabolites
                        metabs = row['Metabolites']
                        stoich = row['stoichiometry']
                        species = metabs.split()
                        coeff = stoich.split()
                        coefficients = []
                        for number in coeff:
                            number = float(number)
                            coefficients.append(number)

                        # Creating the reaction
                        ID = row['Reaction_ID']
                        reaction = Reaction(ID)
                        reaction.name = row['Reaction_name']
                        reaction.subsystem = row['(Reaction_subsystem)']
                        reaction.lower_bound = float(row['Reaction_lb'])
                        reaction.upper_bound = float(row['Reaction_ub'])

                        # Adding metabolites to reaction
                        for n in range(len(species)):
                            m = universal.metabolites.get_by_id(species[n])
                            reaction.add_metabolites({m: coefficients[n]})

                        reaction.gene_reaction_rule = row['Reaction_gene_ (ID_ and/or)']
                        model.add_reactions(reaction)
                        print('\nThe reaction {} has been added to E. coli model'.format(ID))


                    elif row['Reaction_ID'] in in_uni:
                        print(
                            '''\nThe transport reaction for {} contains 
                            metabolites absent in the universal model. The 
                            transport reaction will be added from the input file.'''.format(
                                metab))
                        # The names have to be separated by a space. 
                        # If a metab does not have formula or 
                        # charge the symbol / or 0 respectively 
                        # must be used to indicate it is ot known
                        new_metabs_name = row['metabolite_name'] 
                        new_metab_ID = row['Metabolite_ID']
                        new_metabs_comp = row['Metabolite_compartment']
                        
                        nm_names = new_metabs_name.split()
                        nm_IDs = new_metab_ID.split()
                        nm_comp = new_metabs_comp.split()
                        
                        #print('Name = ', nm_names, 'ID= ', nm_IDs, 'compartment = ', nm_comp)

                        # Reading the metabolites involved in the reaction
                        metabs = row['Metabolites']
                        stoich = row['stoichiometry']
                        species = metabs.split()
                        coeff = stoich.split()
                        
                        coefficients =[]
                        for number in coeff:
                            number = float(number)
                            coefficients.append(number)
                        #print(coefficients)
                            
                        # Creating the reaction
                        ID = row['Reaction_ID']
                        reaction = Reaction(ID)
                        reaction.name = row['Reaction_name']
                        reaction.subsystem = row['(Reaction_subsystem)']
                        reaction.upper_bound = float(row['Reaction_ub'])
                        reaction.lower_bound = float(row['Reaction_lb'])
                        
                        
                        # Create metabolites
                        for n in range(len(nm_names)):
                            metabolite = nm_IDs[n]
                            metabolite = Metabolite(
                                str(nm_IDs[n]), name=str(nm_names[n]), compartment=str(nm_comp[n]))
                            universal.add_metabolites([Metabolite(nm_IDs[n])])
                        
                        
                        # Adding metabolites  already in uni to new reaction
                        for n in range(len(species)):
                            reaction.add_metabolites({model.metabolites.get_by_id(species[n]): coefficients[n]})
                        model.add_reactions(reaction)
                        reaction.gene_reaction_rule = row['Reaction_gene_ (ID_ and/or)']
                        print('\nThe reaction {} has been added to the reference model'.format(ID))
                    else:
                        print('''\nThere is no tranport reaction for {} in 
                        the input file. A transported should be added to the 
                        file before running the analysis'''.format(metab))

def add_reactions_new_metab(input_file, ID, universal, model):
    """
    Add reaction from input file to model. 
    
    First it checks if the reactions have any metab not in universal.
    If so it adds those metabs before adding the desired reaction(s) 
    identified with the ID which must correspond to the one in the 
    column Reaction_ID.
    -------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    ID--str ID of the reaction as found with get_reactions function
    universal--cobra.Model universal model
    model--cobra.Model reference model
    """
    #TODO: delete since unused
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        metabolites = get_metabolites(input_file)
        for row in reader:
            if ID == row['Reaction_ID']:
                # Check if any reaction has metabolite to be added to the universal
                #in_uni = check_if_metabs_in_universal(input_file)

                #if row['Reaction_ID'] in in_uni:  # ?

                    # The names have to be separated by a space. If a metab does not have formula or
                    # charge the symbol / or 0 respectively must be used to indicate it is ot known
                    new_metabs_name = row['metabolite_name']
                    new_metab_ID = row['Metabolite_ID']
                    new_metabs_comp = row['Metabolite_compartment']

                    nm_names = new_metabs_name.split()
                    nm_IDs = new_metab_ID.split()
                    nm_comp = new_metabs_comp.split()

                    # print('Name = ', nm_names, 'ID= ', nm_IDs, 'compartment = ', nm_comp)

                    # Reading the metabolites involved in the reaction
                    metabs = row['Metabolites']
                    stoich = row['stoichiometry']
                    species = metabs.split()
                    coeff = stoich.split()

                    coefficients = []
                    for number in coeff:
                        number = float(number)
                        coefficients.append(number)
                    # print(coefficients)

                    # Creating the reaction
                    ID = row['Reaction_ID']
                    reaction = Reaction(ID)
                    reaction.name = row['Reaction_name']
                    reaction.subsystem = row['(Reaction_subsystem)']
                    reaction.upper_bound = float(row['Reaction_ub'])
                    reaction.lower_bound = float(row['Reaction_lb'])
                    reaction.notes['KEGG Reaction']=row['kegg_id']
                    reaction.notes['EC Number']=row['ec_number']

                    # Create metabolites
                    for n in range(len(nm_names)):
                        metabolite = nm_IDs[n]
                        metabolite = Metabolite(str(nm_IDs[n]), name=str(nm_names[n]), compartment=str(nm_comp[n]))
                        universal.add_metabolites([Metabolite(nm_IDs[n])])

                    # print(list(enumerate(species)), list(enumerate(nm_IDs)))

                    # for index, m in enumerate(nm_IDs):
                    #       reaction.add_metabolites({Metabolite(m): coefficients[index]})

                    # reaction.add_metabolites({Metabolite(str(nm_IDs[i])): str(coefficients[index])})

                    # for ind, metabol in enumerate(species):
                    #   reaction.add_metabolites({Metabolite(i[for i in nm_IDs]):coefficients[ind]})

                    # Adding metabolites  already in uni to new reaction
                    for n in range(len(species)):
                        reaction.add_metabolites({universal.metabolites.get_by_id(species[n]): coefficients[n]})
                    
                    universal.add_reactions(reaction)
                    reaction.gene_reaction_rule = row['Reaction_gene_ (ID_ and/or)']
                    print('\nThe reaction {} has been added to the universal model'.format(ID))

def r_for_uni(input_file, universal, model):
    """
    Checks if the reactions in the input file are present in the reference model. 

    If the reference model does not have them, it checks if they are present 
    in the universal. if they are in the universal and if they are not there 
    the fuction adds them to the universal.
    -------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    universal--cobra.Model universal model
    model--cobra.Model reference model

    Retrun:
    universal: cobra.Model universal model with the added reactions.
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')

        reactions_ids = get_reactions(input_file)
        react_to_be_added_to_uni = []
        #not_in_uni = check_if_metabs_in_universal(input_file)
        for reaction in reactions_ids:
            if reaction in [i.id for i in model.reactions]:
                continue
                # print("\nThis reaction ({}) is already in the E.coli model".format(reaction))
            elif reaction in [i.id for i in universal.reactions]:
                continue
                # print("\nThis reaction ({}) is already in the universal model".format(reaction))
            elif "EX_" not in reaction:
                react_to_be_added_to_uni.append(reaction)
        if react_to_be_added_to_uni != '':
            print('\nThe following reactions will be added to the universal model: {}'.format(react_to_be_added_to_uni))
            for row in reader:
                for ID in react_to_be_added_to_uni:
                    if ID in row['Reaction_ID']:
                        r = cobra.Reaction(ID)
                        universal.add_reaction(r)
                        equation = row['Reaction_equation']
                        r.build_reaction_from_string(equation, fwd_arrow='-->')
                        r.name = row['Reaction_name']
                        r.lower_bound = int(row['Reaction_lb'])
                        r.upper_bound = int(row['Reaction_ub'])
                        r.gene_reaction_rule = row['Reaction_gene_ (ID_ and/or)']
                        r.notes['KEGG Reaction']=row['kegg_id']
                        r.notes['EC Number']=row['ec_number']
                        # print('\nThe following reactions have be added to the universal model: {}'.format(ID))

def set_bounds_ex(input_file, model, metab):
    """
    Set bounds of the exchange reactions of the selected compounds

    If the compounds have to be consumed the bounds of its exchange 
    reaction it set to (-1000, 0).While it it has to be produced they are
    set to (0, 1000)
    -------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model
    metab--str BiGG ID of the metabolites indicated in the input file
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        metabolites = get_metabolites(input_file)
        if metab in metabolites:
            for row in reader:
                    if (row['consumption'] in ('yes', 'Yes', 'Y', 'y', 'YES') and 
                            row['production'] in ('No', 'no', 'NO', 'N', 'n') and 
                            metab in row['metabolite_BiGG_ID']):
                        exchange = model.reactions.get_by_id('EX_' + metab + '_e')
                        print('Exchange {}: '.format(metab), exchange.reaction, 'Old bounds: ', exchange.bounds)
                        exchange.lower_bound = -1000
                        exchange.upper_bound = 0
                        print('Exchange {}: '.format(metab), exchange.reaction, 'New bounds: ', exchange.bounds)
                    elif (row['consumption'] not in ('yes', 'Yes', 'Y', 'y', 'YES') and 
                            row['consumption'] != '' and 
                            row['production'] in ('No', 'no', 'NO', 'N', 'n') and 
                            metab in row['metabolite_BiGG_ID']):
                        exchange = model.reactions.get_by_id('EX_' + metab + '_e')
                        print('Exchange {}: '.format(metab), exchange.reaction, 'Old bounds: ', exchange.bounds)
                        exchange.lower_bound = -1000
                        exchange.upper_bound = 0
                        print('Exchange {}: '.format(metab), exchange.reaction, 'New bounds: ', exchange.bounds)
                    elif (row['production'] in ('yes', 'Yes', 'Y', 'y', 'YES') and 
                            row['consumption'] in ('No', 'no', 'NO', 'N', 'n') and 
                            metab in row['metabolite_BiGG_ID']):
                        exchange = model.reactions.get_by_id('EX_' + metab + '_e')
                        print('Exchange {}: '.format(metab), exchange.reaction, 'Old bounds: ', exchange.bounds)
                        exchange.lower_bound = 0
                        exchange.upper_bound = 1000
                        print('Exchange {}: '.format(metab), exchange.reaction, 'New bounds: ', exchange.bounds)
                    elif (row['production'] not in ('yes', 'Yes', 'Y', 'y', 'YES') and 
                            row['production'] != '' and 
                            row['consumption'] in ('No', 'no', 'NO', 'N', 'n') and 
                            metab in row['metabolite_BiGG_ID']):
                        exchange = model.reactions.get_by_id('EX_' + metab + '_e')
                        print('Exchange {}: '.format(metab), exchange.reaction, 'Old bounds: ', exchange.bounds)
                        exchange.lower_bound = 0
                        exchange.upper_bound = 1000
                        print('Exchange {}: '.format(metab), exchange.reaction, 'New bounds: ', exchange.bounds)
                    #else:
                        #continue

def parser(input_file, universal, model):
    """
    Initialize the reference and universal models for the analysis.

    1) Check if the metabolites are in the medium
    2) Checks if transport reactions are already in the model 
    3) Checks if all the reactions in input are in the model, 
    else checks the universal and if absent it adds them 
    -------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    universal--cobra.Model universal model
    model--cobra.Model reference model

    Retrun:
    model: cobra.Model reference model with the right settings.
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        #fieldnames = ['ï»¿metabolite_BiGG_ID', 'C-source?_y/n', 'N-source?_y/n', 'Reaction_ID', 'Reaction_name',
                      #'(Reaction_subsystem)', 'Reaction_equation', 'Reaction_lb', 'Reaction_ub',
                      #'Reaction_gene_ (ID_ and/or)', 'Metabolites', 'stoichiometry', 'metabolite_name',
                      #'Metabolite_compartment', 'Metabolite_formula', 'Metabolite_charge']
        reader = csv.DictReader(csvfile, dialect='excel')
        metabolites = get_metabolites(input_file)  # Function
        # check if the metabolites is in the medium and if not add exchange
        
            #elif row['production'] != "":
                
                #if ('EX_' + i + '_e') in model.reactions:
                #    print('{} is in the medium'.format(i), '\n')
                #elif ('EX_' + i + '_e') not in model.reactions:
                #    print('{} is not in the medium'.format(i), '\n')
                #    if ('EX_' + i + '_e') in universal.reactions:
                #        exchange = universal.reactions.get_by_id('EX_' + i + '_e')
                #        exchange.lower_bounds = 0
                #        exchange.upper_bounds = 1000
                #        model.add_reaction(exchange)
                #    else:
                #        extracellular = Metabolite(i + '_e')
                #        extracellular = Metabolite(i + '_e', name='extracellular' + i, compartment='_e')
                #        exch = Reaction('EX_' + i + '_e')
                #        exch.name = 'Exchange {}'.format(extracellular)
                #        exch.lower_bound = 0
                #        exch.upper = 1000
                #        exch.add_metabolites({extracellular: -1})
                #        model.add_reaction(exch)
                #        print('\nTher reaction {} has been added to E.coli model, hence {} is now in the medium'.format(
                #            'EX_' + i + '_e', i + '_e'))

        # Check for transport reactions in the reference model 
        for metab in metabolites:
            solution = []
            for r in model.reactions:
                if (' '+metab + "_e") in r.reaction and (' '+metab + '_c') in r.reaction:
                    uptake = model.reactions.get_by_id(r.id)
                    if uptake.get_coefficient((metab + '_c')) > 0:
                        solution.append((r.id, r.reaction))
                        print(
                            '\nFor {} there is already a transport reaction allowing the uptake from the extracellular space: '.format(
                                metab),
                            '\nReaction ID: ', r.id,
                            '\nReaction equation: ', r.reaction)
                elif (' '+metab + "_p") in r.reaction and (' '+metab + '_c') in r.reaction:
                    logging.debug(r.id)
                    uptake = model.reactions.get_by_id(r.id)
                    if uptake.get_coefficient(metab + '_c') > 0:
                        solution.append((r.id, r.reaction))
                        print(
                            '\nFor {} there is already a transport reaction allowing the uptake from the periplasmic space: '.format(
                                metab),
                            '\nReaction ID: ', r.id,
                            '\nReaction equation: ', r.reaction)
                else:
                    continue
            if len(solution) == 0:  # if reactions not in reference model checks if they are in universal.
                print("For {} there isn't any uptake trasnsporter in the reference model".format(metab))  # FIX/IMPROVE!
                transporter = 0
                for i in universal.reactions:
                    if (metab + "_e") in i.reaction and (metab + '_c') in i.reaction:
                        uptake = universal.reactions.get_by_id(i.id)
                        if uptake.get_coefficient((metab + '_c')) > 0:
                            print(
                                '\nFor {} there is a transport reaction in the universal model for the uptake from the extracellular space: '.format(
                                    metab),
                                '\nReaction ID: ', i.id,
                                '\nReaction equation: ', i.reaction)
                            model.add_reaction(uptake)
                            transporter+=1
                    elif (metab + "_p") in i.reaction and (metab + '_c') in i.reaction:
                        uptake = universal.reactions.get_by_id(i.id)
                        if uptake.get_coefficient((metab + '_c')) > 0:
                            print(
                                '''\nFor {} there is a transport reaction 
                                    in the universal model for the uptake 
                                    from the periplasm: '''.format(metab),
                                '\nReaction ID: ', i.id,
                                '\nReaction equation: ', i.reaction, '\n')
                            model.add_reaction(uptake)
                            transporter+=1

                # elif [(metab+'e') in row['Reaction_equation'] and (metab+'c') in row['Reaction_equation']] or [(metab+'p') in row['Reaction_equation'] and (metab+'c') in row['Reaction_equation']]:
                #   add_transport_with_metabs_already_in_universal(input_file, metab)

                if transporter == 0:  # if they are not in the universal it checks if they are in the input file and if so add them
                    for row in reader:
                        if (metab+'_e' in row['Reaction_equation'] and
                            metab+'_p' in row['Reaction_equation']) or(
                            metab+'_e' in row['Reaction_equation'] and
                            metab+'_c' in row['Reaction_equation']) or(
                            metab+'_p' in row['Reaction_equation'] and
                            metab+'_c' in row['Reaction_equation']):
                                r = cobra.Reaction(row['Reaction_ID'])
                                model.add_reaction(r) #the transporter should be added default to the model
                                print("""The trasporter has been added to the 
                                reference model from the input file""")
                                equation = row['Reaction_equation']
                                r.build_reaction_from_string(equation, fwd_arrow='-->')
                                r.name = row['Reaction_name']
                                r.lower_bound = int(row['Reaction_lb'])
                                r.upper_bound = int(row['Reaction_ub'])
                                r.gene_reaction_rule = row['Reaction_gene_ (ID_ and/or)']
                                r.notes['KEGG Reaction']=row['kegg_id']
                                r.notes['EC Number']=row['ec_number']
                    #add_transport_with_metabs_already_in_universal(input_file, metab, universal, model)  # Function
                    # Signal somehow that it should be added to the input file?
                    #
        r_for_uni(input_file, universal, model)
        
    return model

#TODO: use cobra function to build the reactions from string