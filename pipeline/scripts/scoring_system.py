#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for score and compare the difference metabolic engineering strategies

Contains functions to evaluate the engineering strategies on 
a set of criteria: consumption rate, production rate, number of 
reactions Knock ins, thermodynamic driving force and pathway
length. The scores are calculated with an additive scoring function
using normalized criteria'scores and user-selected weights.  
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"

import logging
import sys

import cobra
import csv
import numpy
from matplotlib import pyplot as plt

from pipeline.scripts.analysis import *

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def read_weights(input_file):
    """
    Read criteria's weights from input file

    Read the weight that the user has assigned to each criterion.
    The user should type the criteria under 'criteria' column in the
    input file. Some name should be added default:
    - interventions
    - mdf
    - path_length
    While consumption and production should be written with the 
    following format: consumption_+metabolite ID 
    (or production_+metabolite ID).
    'consumption_' or 'production_' must be followed by the metabolite 
    ID (without compartment) as written below 'metabolite_BiGG_ID' 
    column.

    The weights will be normed and their sum equal 1. 
    -----------------------------------------------------------------
    Argument:
    input_file--str, input file in .csv format dictionary like
    
    Return:
    dict_weights: dict with the criteria as keys and their normed
        weights as values
    """
    to_consume = consumption_metabs(input_file)
    to_produce = get_production_objectives(input_file)
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        dict_weights ={}
        for row in reader:
            for m in to_consume:
                if row['criteria'] == 'consumption_' + m:
                    if row['weight'] != '':
                        dict_weights[m] = row['normed_weight']
                    #else:
                    #    dict_weights[m] = 1
            for p in to_produce:
                if row['criteria'] == 'production_' + p:
                    if row['weight'] != '':
                        dict_weights[p] = row['normed_weight']
                    #else:
                    #    dict_weights[p] = 1
                    
            if row['criteria'] == 'interventions':
                if row['normed_weight'] != '':
                    dict_weights['interventions'] = row['normed_weight']
                #else:
                #    dict_weights['interventions'] = 1
                
            elif row['criteria'] == 'mdf':
                if row['normed_weight'] != '':
                    dict_weights['MDF'] = row['normed_weight']
                #else:
                #    dict_weights['MDF'] = 1

            elif row['criteria'] == 'path_length':
                if row['normed_weight'] != '':
                    dict_weights['path_length'] = row['normed_weight']
                #else:
                #    dict_weights['path_length'] = 1

    return dict_weights

def get_max_and_min_values(input_file, output_consumption, output_con_and_prod):
    """
    Collect values of all models for each criterion. 

    Read the output files and collects the values of each critera
    from each model variant. 
    -----------------------------------------------------------------
    Argument:
    input_file--str, input file in .csv format dictionary like
    output_consumption--dict retured by analysis_gf_sol function
    output_con_and_prod--dict retured by cons_prod_dict function
        
    Return:
    consumption_values: dict with metabolites as key and list of 
        consumption rates as values 
    production_values: dict with metabolites as key and list of 
        production rates as values
    number_interventions: list of number of reaction knock-ins of 
        each model variant
    mdf: list of MDF values of each model variant
    p_length: list of number of reactions involved in the pathway 
        of each model variant
    """
    # modify 
    # instead of output_consumption it is the ele with index [6] in output_con_and_prod
    to_consume = consumption_metabs(input_file)
    consumption_values = {}
    for compound in to_consume:
        values_metab = []
        for n in range(1, len(output_consumption)+1): 
            #print('---------'+str(n)+'---------')                                    
        
            
            for key, value in output_con_and_prod['production_'+str(n)].items():
                # if condition that GF production not needed 
                if '_Run_' not in key:
                    dict_cons = output_con_and_prod['consumption_'+str(n)][1][3]
                    for k, v in dict_cons.items():
                        if compound in k:
                            if v not in values_metab:
                                values_metab.append(v)

                # GF needed --> then read the [6]
                else: 
                    consumption_rate = value[6]
                    for k, v in consumption_rate.items():
                        if compound in k:
                            if v not in values_metab:
                                    values_metab.append(v)
        consumption_values[compound] = values_metab

    
    mdf = []
    p_length = []
    for n in range(1, len(output_consumption)+1):
        #print('---------'+str(n)+'---------')
        for i, values in output_con_and_prod['production_'+str(n)].items():
            logging.debug(i, values)
            if '_Run_' not in i:
                # read the thrmodynamic values
                thermo = values['thermodynamic']
                logging.debug(thermo)
                if thermo['mdf']!=None:
                    mdf.append(float(thermo['mdf']))
                else:
                    mdf.append(0)

                if thermo['pathway_length']!=None:
                    p_length.append(float(thermo['pathway_length']))
                else:
                    p_length.append(0)
            else:
                thermo = values[5]
                if thermo['mdf']!=None:
                    mdf.append(float(thermo['mdf']))
                else:
                    mdf.append(0)

                if thermo['pathway_length']!=None:
                    p_length.append(float(thermo['pathway_length']))
                else:
                    p_length.append(0)
    


    to_produce = get_production_objectives(input_file)
    production_values = {}
    for metab in to_produce:
        pvaluesmetab = []
        for n in range(1, len(output_consumption)+1):
            #print('---------'+str(n)+'---------')
            for key, value in output_con_and_prod['production_'+str(n)].items():
                # if condition that GF production not needed 
                if '_Run_' not in key:
                    for k, val in value.items():
                        if metab in k:
                            print(val)
                            if val not in pvaluesmetab:
                                pvaluesmetab.append(val)
                                
                else:
                    prod = value[3]
                    for k, val in prod.items():
                        if metab in k:
                            if val not in pvaluesmetab:
                                pvaluesmetab.append(val)

        production_values[metab] = pvaluesmetab
    
        
    # number of interventions (knock-ins)
    number_interventions = []
    for n in range(1, len(output_consumption)+1):
            #print('---------'+str(n)+'---------')
            for key, values in output_con_and_prod['production_'+str(n)].items():
                if '_Run_' in key:
                    dict_KI = values[4]
                    for entry in list(dict_KI):
                        if dict_KI[entry] == 0: # eliminate the reactions with a zero flux
                            dict_KI.pop(entry)
                    length = len(dict_KI)
                else:
                    list_KI_cons = output_con_and_prod['consumption_'+str(n)][1][0]
                    length = len(list_KI_cons)
                number_interventions.append(length) 
                #for m in range(1, len(output_con_and_prod['production_'+str(n)])+1):
                #    dict_KI = output_con_and_prod['production_'+str(n)]['lac__L_Run_'+str(m)][4]
                #    for entry in list(dict_KI):           
    return consumption_values, production_values, number_interventions, mdf, p_length


def generate_dict_ranges_per_criteria(ranges):
    """
    Generate tuple range with max and min values of each criterion

    Given the list of values for each criteria and each variant, it
    identifies the maximum and minumum values that will be used for 
    the normalization.
    -----------------------------------------------------------------
    Argument:
    ranges--dictionaries and list returned by get_max_and_min_values
        function
        
    Return:
    values_range: dict with the metaolites or the criteria as keys 
        and the tuples (min, max) values as values 
        e.g. {metabA:(min_a, max_a), metabB:(min_b, max_b), 
            'Interventions': (min_c, max_c)
            }
    """
    values_range = {}
    
    for key, value in ranges[0].items():
        logging.debug('cons values')
        logging.debug(value)
        if max(value) == min(value):
            values_range[key] = (min(value), 0) #minimized
            print('For {} the normalization should be done with {}'.format(key, values_range[key]))
            
        else:
            print('For {} the normalization has to be done using {} as min and {} as max value'.format(key, min(value), max(value)))
            values_range[key] = (min(value), max(value))
    
    for key, value in ranges[1].items():
        logging.debug('prod values')
        logging.debug(value)
        if max(value) == min(value):
            values_range[key] = (0, max(value)) #maximized
            print('For {} the normalization should be done with {}'.format(key, values_range[key]))
             
        else:
            print('For {} the normalization has to be done using {} as min and {} as max value'.format(key, min(value), max(value)))
            values_range[key] = (min(value), max(value))
            
    if max(ranges[2]) == min (ranges[2]):
        values_range['Interventions'] = (min(ranges[2]), 0) #minimized
        print('No normalization is needed for the number of heterologous reactions: {}'.format(values_range['Interventions']))
        
    else:
        values_range['Interventions'] = (min(ranges[2]), max(ranges[2]))

    # add condition for if GF not needed
    # mdf not calculated if GF not needed - fix that.

    if max(ranges[3]) == min (ranges[3]):
        values_range['MDF'] = None #maximized
        print('No normalization is needed for the mdf values: {}'.format(values_range['MDF']))
        
    else:
        values_range['MDF'] = (min(ranges[3]), max(ranges[3]))

    if max(ranges[4]) == min (ranges[4]):
        values_range['path_length'] = None #minimized
        print('No normalization is needed for pathway length: {}'.format(values_range['path_length']))
    else:
        values_range['path_length'] = (min(ranges[4]), max(ranges[4]))
        
    return values_range


def generate_scores(input_file, output_consumption, output_con_and_prod):
    """
    Calculate scores with additive scoring function. 
    
    Calculate the scores and writes the partof_summary_output.csv file
    inside pipeline/outputs folder. This file will be used to
    initiate the summary_output.csv file with the scores of each 
    model variant.
    Criteria that are preferred when have higher values (i.e. MDF and 
    production rate) are scored with the following formula:
        (Value - Value_min) / (Value_max - Value_min)

    Whereas, to score criteria that should be minimised (i.e. number of 
    heterologous, consumption rate, pathway length) this other formula 
    is used:
        (Value_max - Value) / (Value_max - Value_min)

    The final score for each model is obtained by calculating the 
    weigthed sum for the scores for each criteria[1]. 

    [1] Schneider, P., & Klamt, S. (2019). Characterizing and 
        ranking computed metabolic engineering strategies. 
        Bioinformatics, 35(17), 3063–3072. 
        https://doi.org/10.1093/bioinformatics/bty1065
    -----------------------------------------------------------------
    Arguments:
    input_file--str, input file in .csv format dictionary like
    output_consumption--dict retured by analysis_gf_sol function
    output_con_and_prod--dict retured by cons_prod_dict function
    
    Returns:
    scores_list: list with the scores accessible with index 0 in the 
        called function
    score_dict: dict with the scores and the associated models
        accessible with index 1 in the called function
    """
    #TODO: add min and max reference values for MDF and path_length

    dict_w = read_weights(input_file)
    criteria=[]
    weights=[]
    for crit, weight in dict_w.items():
        criteria.append(crit)
        weights.append(weight)

    logging.debug(dict_w)
    to_consume = consumption_metabs(input_file)
    to_produce = get_production_objectives(input_file)
    ranges = get_max_and_min_values(input_file, output_consumption, output_con_and_prod)
    logging.debug(ranges)
    ranges_per_criteria = generate_dict_ranges_per_criteria(ranges)
    logging.debug(ranges_per_criteria)

    with open('pipeline/outputs/partof_summary_output.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["model_name", 
                        "criteria", 
                        "weights", 
                        "normalised_scores", 
                        "tot_metabolic", 
                        "percentage_Biobrick",
                        "allowed_KOs",
                        "reactions_KOs",
                        "objective_value"])
        writer.writerow(["#identifier of the model: the first number indicate the variant consumption's reaction search, this is followed by the identifier of the target, finally, the last number indicates the variant for production's reaction search", 
                        "#criteria used to compare the model variants between each other. Default ones are fluxk thorught the exchange reaction of the substrate (minimized), production rate fo the target (maximized), number of reactions knock ins (minimized), MDF value (maximised), length of the pathway (minimized)", 
                        "#user selected weights for each criteria as assigned undr the column 'weight' in the input file", 
                        "#normalized scores for each criteria. Each normalized score ranges from 0 to 1, with 1 indicating the best among the model variants", 
                        "#final score of each model. Derived by multiplying the normed scores by the weight of theri criteria. The final scores ranges from 0 to 1. Scores closer to 1 indicate the most favourable engineering strategies according to the computational metabolic analysis", 
                        "#percentage indicating the number of biobricks found per number of reactions knock ins suggested",
                        "#number of reaction knock outs allowed per optknock analysis round",
                        "#BiGG reactions identifiers for the knock outs found with optknock analysis.",
                        "#objective value resulting from running FBA on the model when the knock-outs are applied. The objective of this optimization is the production rate of the target"])
        #for key, value in ranges[1].items():
        #    if max(value) == min(value):
        #        print('For {} the normalization should be done with -1000 and 1000'.format(key))
        #    else:
        #        print('For {} the normalization has to be done using {} as min and {} as max value'.format(key, min(value), max(value)))

        scores_list = []
        score_dict = {}

        #for key, value in ranges[0].items():
        #    if max(value) == min(value):
        #        print('For {} the normalization should be done with -1000 and 1000'.format(key))
        #    else:
        #        print('For {} the normalization has to be done using {} as min and {} as max value'.format(key, min(value), max(value)))
        list_single_scores = []
        for n in range(1, len(output_consumption)+1):
                print('---------'+str(n)+'---------')
                for key, values in output_con_and_prod['production_'+str(n)].items():
                #for m in range(1, len(output_con_and_prod['production_'+str(n)])+1):
                    if '_Run_' in key:
                        dict_KI = values[4]
                        
                        for entry in list(dict_KI):
                            if dict_KI[entry] == 0:
                                dict_KI.pop(entry)
                        number_of_interventions = len(dict_KI)


                    else:
                        list_KI_cons = output_con_and_prod['consumption_'+str(n)][1][0]
                        number_of_interventions = len(list_KI_cons)
                        
                    
                    if '_Run_' in key:
                            dict_production = values[3].values()
                            for i in dict_production:
                                production_val = i
                            logging.debug(production_val)
                    else:
                        target = to_produce[0]
                        production_val = values['EX_'+target+'_e flux']
                    

                    logging.debug('consumption')
                    scores_consumption={}
                    if '_Run_' not in key:
                        dict_consumption = output_con_and_prod['consumption_'+str(n)][1][3]
                        for c in to_consume:
                                range_tuple = ranges_per_criteria[c]
                                print(c, range_tuple)
                                value = dict_consumption['EX_'+c+'_e']
                                score = (float(range_tuple[1])-float(value))/(float(range_tuple[1])-float(range_tuple[0]))
                                scores_consumption[c]=score
                                list_single_scores.append(score)
                                logging.debug(score)
                    else:
                        dict_consumption = values[6]
                        for c in to_consume:
                                range_tuple = ranges_per_criteria[c]
                                print(c, range_tuple)
                                value = dict_consumption['EX_'+c+'_e']
                                score = (float(range_tuple[1])-float(value))/(float(range_tuple[1])-float(range_tuple[0]))
                                scores_consumption[c]=score
                                list_single_scores.append(score)
                                logging.debug(score)
                        
                    logging.debug('production')
                    scores_production={}
                    for p in to_produce:
                            range_tuple = ranges_per_criteria[p]
                            print(p, range_tuple)
                            value = production_val
                            score = (value-(range_tuple[0]))/(range_tuple[1]-range_tuple[0])
                            scores_production[p] = score
                            list_single_scores.append(score) 
                            logging.debug(score)      
                        
                    logging.debug('interventions')
                    if len(values) == 9:
                            number_of_interventions = number_of_interventions + values[8]
                        #logging.debug('Number of Interventions', ranges_per_criteria['Interventions'])
                    if ranges_per_criteria['Interventions'] == None:
                            Sh = 1
                            print('score_c: ', scores_consumption, 'score_p: ', scores_production, 'score_KI: ', Sh)
                            score_total = float(Sh)
                    else:
                            logging.debug(ranges_per_criteria['Interventions'][1])
                            Sh = (int(ranges_per_criteria['Interventions'][1]) - number_of_interventions)/(int(ranges_per_criteria['Interventions'][1])-int(ranges_per_criteria['Interventions'][0]))
                            print('score_c: ', scores_consumption, 'score_p: ', scores_production, 'score_KI: ', Sh)
                            Wh = dict_w['interventions']
                            logging.debug(Wh)
                            score_total = float(Sh)*float(Wh)
                    list_single_scores.append(Sh)

                    logging.debug(score_total)

                    for x in scores_consumption.keys():

                            score_total += float(scores_consumption[x])*float(dict_w[x])

                    for y in scores_production.keys():
                            score_total += float(scores_production[y])*float(dict_w[y])

                    if '_Run_' in key:
                        #TODO: in function for final dict add considering mdf and path lenght even if no KI are needed for production
                        #---- INSERT VALUE MDF -----
                        mdf = values[5]['mdf']
                        logging.debug('mdf')
                        if mdf != None:
                            range_val = ranges_per_criteria['MDF']
                            Smdf = (float(mdf) -(range_val[0]))/(range_val[1]-(range_val[0])) 
                            list_single_scores.append(Smdf)
                            score_total += (float(Smdf)*float(dict_w['MDF']))
                            logging.debug(Smdf)
                            logging.debug(score_total)
                        else:
                            mdf = 0
                            list_single_scores.append(mdf)
                            score_total += (float(mdf)*float(dict_w['MDF']))
                            logging.debug(score_total)

                        #-----INSERT PATHWAY LENGTH ----
                        path_length = values[5]['pathway_length']
                        logging.debug('pl')
                        if path_length != None:
                            range_val = ranges_per_criteria['path_length']
                            Spath_length = (range_val[1]-path_length)/(range_val[1]-(range_val[0])) 
                            list_single_scores.append(Spath_length)
                            score_total += float(Spath_length)*float(dict_w['path_length'])
                            logging.debug(Spath_length)
                            logging.debug(score_total)
                        else:
                            path_length = 0
                            list_single_scores.append(path_length)
                            score_total += float(path_length)*float(dict_w['path_length'])
                            logging.debug(score_total)
                    else:
                        mdf = values['thermodynamic']['mdf']
                        logging.debug('mdf')
                        if mdf != None:
                            range_val = ranges_per_criteria['MDF']
                            Smdf = (float(mdf) -(range_val[0]))/(range_val[1]-(range_val[0])) 
                            list_single_scores.append(Smdf)
                            score_total += (float(Smdf)*float(dict_w['MDF']))
                            logging.debug(Smdf)
                            logging.debug(score_total)
                        else:
                            mdf = 0
                            list_single_scores.append(mdf)
                            score_total += (float(mdf)*float(dict_w['MDF']))
                            logging.debug(score_total)

                        #-----INSERT PATHWAY LENGTH ----
                        path_length = values['thermodynamic']['pathway_length']
                        logging.debug('pl')
                        if path_length != None:
                            range_val = ranges_per_criteria['path_length']
                            Spath_length = (range_val[1]-path_length)/(range_val[1]-(range_val[0])) 
                            list_single_scores.append(Spath_length)
                            score_total += float(Spath_length)*float(dict_w['path_length'])
                            logging.debug(Spath_length)
                            logging.debug(score_total)
                        else:
                            path_length = 0
                            list_single_scores.append(path_length)
                            score_total += float(path_length)*float(dict_w['path_length'])
                            logging.debug(score_total)
                    scores_list.append(score_total)
                    #print('\n', score_total)
                    tot_metab_score = [score_total]
                    score_dict[str(n)+'_'+key] = score_total
                    for i in range(len(criteria)):
                            if (i+1) <= len(tot_metab_score):
                                writer.writerow([(str(n)+'_'+key), criteria[i], weights[i], list_single_scores[i], tot_metab_score[i]])
                            elif (i+1) > len(tot_metab_score):
                                writer.writerow(['', criteria[i], weights[i], list_single_scores[i]])
        
        writer.writerow(['-----------------------------------------------------'])
        writer.writerow(['rank_position', 'model_name', 'final_metabolic_score'])
        writer.writerow(['#position on the model variants once ranked by the final scores in decreasing order', 
                        "#identifier of the model: the first number indicate the variant consumption's reaction search, this is followed by the identifier of the target, finally, the last number indicates the variant for production's reaction search", 
                        "#final score of each model. Derived by multiplying the normed scores by the weight of theri criteria. The final scores ranges from 0 to 1. Scores closer to 1 indicate the most favourable engineering strategies according to the computational metabolic analysis"])
                        #summary_output(n, dict_w, (str(n)+'_'+key), list_single_scores, tot_metab_score)
    return scores_list, score_dict

def scores_evaluations(input_file, output_consumption, output_con_and_prod):
    """
    Rank the different engineering strategies

    This function parses the output of generate_scores function and 
    ranks the model variants representing different engineering 
    strategies by total scores in decreasing order. The model at 
    position 1 is the one considered the best based on the 
    "metabolic criteria": consumption and production rate, number of
    knock-ins, MDF,pathway length. 
    -----------------------------------------------------------------
    Arguments:
    input_file--str, input file in .csv format dictionary like
    output_consumption--dict retured by analysis_gf_sol function
    output_con_and_prod--dict retured by cons_prod_dict function

    Return:
    final_output_scores: dict with the ranked variant. The rank
        positons are the keys, while the values are tuples with 
        model name and its finals score.    
    """
    scoring_info = generate_scores(input_file, output_consumption, output_con_and_prod)
    logging.debug(scoring_info)
    scores_evaluation = numpy.array(scoring_info[0])
    sor = (-scores_evaluation).argsort()
    ranksm = numpy.empty_like(sor)
    ranksm[sor] = numpy.arange(len(scores_evaluation))
    final_rank = ranksm+1
    final_rank = ranksm+1
    final_output_scores={}
    keys = scoring_info[1].keys()              
    for n in range(1, len(final_rank)+1):
        for index, rank in enumerate(final_rank):
            if n == rank:
                for i, model in enumerate(keys):
                    if index == i:
                        print("The final position in the ranking for model {} is {} and its final score is {}".format(model, rank, scoring_info[0][i]))
                        final_output_scores[n]= (model, scoring_info[0][i])
                
    
    return final_output_scores


def score_BB_presence(input_file, to_wiki, from_wiki, scores_output):
    """
    Get percentage of matching biobrick per model variant

    Counts the number of knock ins involved in each different
    engineering strategy (model variant) and it counts the number
    of biobrick matchning the reaction function were found by
    querying the wikidata biobrick database.

    It results the percentage of matching biobricks.
    -----------------------------------------------------------------
    Arguments:
    input_file--str, input file in .csv format dictionary like
    to_wiki--list of lists with information on EC number
        and KEGG ID of each heterologous reaction
    from wiki--dict type (nested dictionary) returned from 
        wikidata query function.
    scores_output--dict with the ranked variant. The rank
        positons are the keys, while the values are tuples with 
        model name and its finals score.

    Return:
    BB_scores: list of tuples in which the object with index 0 is
        the name of the model variant, while index 1 indicate 
        the percentage of biobricks   
    """
    # Counts the number of heterologous reactions (Knock Ins)
    n=0
    tot_heterologous = []
    for i in to_wiki:
        n+=1
        logging.debug('---{}---'.format(n))
        n_heterologous_reactions = len(i)
        #for reaction_info in i.values():
            #reaction_info['EC']
        logging.debug(n_heterologous_reactions)
        tot_heterologous.append(n_heterologous_reactions)
    # Counts the number of Biobricks (BB) matching the heterologous reactions
    n=0
    tot_BB = []
    for j in from_wiki:
        n+=1
        logging.debug('---{}---'.format(n))
        n_BB = len(j)
        for key, subdict in j.items():
            if subdict['BB_name'] == None:
                logging.debug('What a shame! no BB for {}'.format(key))
                n_BB = n_BB - 1
        logging.debug(n_BB)
        tot_BB.append(n_BB)

    BB_scores = []
    # Calculate score and add to the scoring dictionary 
    for position, tup in scores_output.items():
        raw_fraction = tot_BB[position-1] / tot_heterologous[position-1]
        #normalized = (raw_fraction - minval) / (maxval - minval)
        #logging.debug(normalized)
        old_score = tup[1]
        logging.debug(old_score)
        BB = (tup[0], raw_fraction) # maybe the raw one is more informative, it inform on the percentage, the other informs on the percentage but is relative to the max and min % among the models. None informs on the number of BB.
        BB_scores.append(BB)
    return BB_scores

def plotting_KI_BB_scores(scores_output, scores_bb):
    """
    Visualization of the metabolic vs biobrick scores

    Uses matplotlib to plot the scores of each model variant
    as obtained by evaluating the "metabolic 
    criteria" (i.e. consumption and production rate, number of
    knock-ins, MDF,pathway length) vs the percentage of 
    matching biobricks found per each variant.

    The metabolic score are on the y asix while the % of 
    biobricks are on the x axis
    This function is used only for analysis run on python 
    interpreter, not via command line calls of the analysis.
    -----------------------------------------------------------------
    Arguments:
    scores_output--dict with the ranked variant. The rank
        positons are the keys, while the values are tuples with 
        model name and its finals score.
    scores_bb--list returned by score_BB_presence function.
    """
    metab_scores = []
    bb_scores = []
    labels=[]

    for position, tup in scores_output.items():
        metab_scores.append(tup[1])
    for i in scores_bb:
        bb_scores.append(i[1])
        full_name = i[0]
        short = '  Model '+full_name[0:1]+'.'+full_name[-1]
        labels.append(short)
    print(metab_scores)
    print(bb_scores)
    fig, ax = plt.subplots()
    ax.scatter(bb_scores, metab_scores)

    for i, name in enumerate(labels):
        ax.annotate(name, (bb_scores[i],metab_scores[i]))
    plt.xlabel("% matching Biobricks")
    plt.ylabel('Metabolic score')

