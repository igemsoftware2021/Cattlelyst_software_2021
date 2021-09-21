
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Optknock calls

Module with functions for running Optknock analysis 
to find reactions knock outs that could couple growth to 
production of a target compound.   
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"


# Import statements

import cobra
import logging
import sys
import csv

import pandas
from cobra import Model, Reaction, Metabolite

from cobra.util.solver import check_solver_status
from cobra.exceptions import (
    OPTLANG_TO_EXCEPTIONS_DICT, Infeasible, OptimizationError, SolverNotFound)
 
from pipeline_package.call_Optknock_robustknock import list_excluded_reactions                 
from pipeline_package.Optknock_robustknock import Run_OptKnock
from pipeline_package.analysis import (remove_rlist, get_biomass_equation, 
                    get_ex_c_source_metab, get_production_objectives)


logging.basicConfig(stream=sys.stderr, level=logging.ERROR)

def add_Optknock_to_analysis(input_file): # new function added by Delielena Poli in date 23/09/2020
    """
    Reads the input file to decide whether doing Optknock analysis or not
    -------------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like

    Return:
    optknock: boolean, True if knock outs should be considered in the
        engineering stratiegy, False if KOs should not be considered.
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        for row in reader:
            if row['Optknock'] in ('yes', 'Yes', 'Y', 'y', 'YES'):
                optknock = True
                break
            elif row['Optknock'] in ('No', 'no', 'NO', 'N', 'n'):
                optknock = False
                break

    return optknock


def run_optknock_analysis(output_consumption, output_con_and_prod, model, universal):
    """
    Find KO to couple production and growth.

    Optknock is used to look for reactions knock outs that could
    allow the coupling of production to growth. 
    -----------------------------------------------------------------
    output_consumption--dict as output from analysis_gf_sol function
    output_con_and_prod--dict retured by cons_prod_dict function
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace

    Return:
    per_cons_var: dict with the solutions of Optknock per 
        variant with optimized consumption. 
    """
    per_cons_var={}
    candidates = list_excluded_reactions(model)
    for n in range(1, len(output_consumption)+1):
        print('---------'+str(n)+'---------')
        result_opt = {}
        for key, values in output_con_and_prod['production_'+str(n)].items():
            if type(values) != dict:
                KI_model = values[7]
                #for exchange, flux in values[3].items():
                #    target = exchange
                #KI_list = []
                #for i in KI.keys():
                #    KI_list.append(i)
                #    r = universal.reactions.get_by_id(i)
                #    model.add_reaction(r)
                #logging.debug(KI_list)
                for exchange, flux in values[3].items():
                    target = exchange
                opt = []
                j=1
                dict_per_key = {}
                while j <= 6:
                    print('Running OptKnock . . .')
                    knock, info_dictionary = Run_OptKnock(KI_model,target, j, KO_cand = candidates, verbose=True)
                    opt += [knock] # it should already print stuff but I want also to recod the info in a dictionary/model
                    dict_per_key[j] = info_dictionary
                    print(dict_per_key)
                    j = j+1
                result_opt[str(n)+'_'+key] = (opt, dict_per_key) #opt object for future production envelope, info dict
                print('res_opt =', result_opt)  
                #remove_rlist(KI_list, model)
            else:
                target = 'EX_'+key+'_e'
                print(target)
                opt = []
                j=1
                dict_per_key = {}
                while j <= 6:
                    print('Running OptKnock . . .')
                    knock, info_dictionary = Run_OptKnock(model,target, j, KO_cand = candidates, verbose=True)
                    opt += [knock] # it should already print stuff but I want also to recod the info in a dictionary/model
                    dict_per_key[j] = info_dictionary
                    print(dict_per_key)
                    j = j+1
                result_opt[str(n)+'_'+key] = (opt, dict_per_key) #opt object for future production envelope, info dict
                print('res_opt =', result_opt)  
                #remove_rlist(KI_list, model)
        per_cons_var[n]=result_opt
    return per_cons_var

def production_env_KO_eval(input_file, model, universal, output_optknock, output_cons_prod):
    """
    Do production envelope analysis on mutants with reaction knockouts

    Given the output of Optknock analysis, it evaluate the solutions
    by calculating the production envelope of the model variant with
    the reactions knock-ins and the variants that also include reaction
    knock-outs (KOs).
    
    The output is a the list of fluxes thought the exchange reaction of the
    target when the biomass reaction is increasingly limited in its upper
    bound. These results are obtained for each max allowed knock-out and for
    the model variant with only reaction addition too. 

    This information is also provided in the form of txt files found in the
    pipeline folder.
    -----------------------------------------------------------------
    input_file--str input file in .csv format dictionary like
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace
    output_optknock--dict with the solutions of Optknock per 
        variant with optimized consumption. Returend by 
        run_optknock_analysis function.
    output_con_and_prod--dict retured by cons_prod_dict function

    Return:
    production_envelope_per_model: dict with the pandas.DataFrame 
        objects of the production envelope of each model variants
        and the correspondant KO mutants. 
    """
    biomass = get_biomass_equation(model) #added by Delielena Poli in date 23/09/2020
    c_source = get_ex_c_source_metab(input_file, model) #added by Delielena Poli in date 23/09/2020
    target = get_production_objectives(input_file) #added by Delielena Poli in date 07/10/2020
    #Production envelope analysis added by Delielena Poli in date 07/10/2020
    production_envelope_per_model={}
    
    Niter = 10
    for key, values in output_optknock.items():
        print('-----{}-----'.format(key))
        for i in values: #model variant identifier
            print('-----{}-----'.format(i))
            #print(output_cons_prod['production_'+i[0]].keys())
            print(i[2:])
            for k in output_cons_prod['production_'+i[0]].keys():
                print("k = ", k)
                
                if i[2:] == k:
                    #production envelope with only KIs
                    #create dictionary for production envelope data

                    prod_env_KI_KO = []
                    p_KI = {}
                    try:
                        growth=[]
                        productMax=[]
                        productMin=[]
                        print('KI model')
                        KI_model = output_cons_prod['production_'+i[0]][i[2:]][7]
                        biomass = get_biomass_equation(KI_model)
                        KI_model.objective=biomass 
                        #model.reactions.get_by_id(source[0]).lower_bound=-1 
                        maxG=KI_model.optimize().objective_value
                        print('gmax', maxG)
                        KI_model.objective='EX_'+target[0]+'_e'
                        for o in range(Niter+1):
                            g=maxG/Niter*o
                            print(o, g)
                            biomass.lower_bound=g
                            maxp = KI_model.optimize().objective_value
                            minp = KI_model.optimize(objective_sense='minimize').objective_value
                            productMax.append(maxp)
                            productMin.append(minp)
                            growth.append(g)
                            print(maxp)
                            print(minp)


                            with open("Prod_envelope_Aero"+i+".txt", 'w') as f:
                                for item in growth:
                                    f.write("%s\n" % item)


                            with open("Prod_envelope_Aero_Max"+i+".txt", 'w') as f:
                                for item in productMax:
                                    f.write("%s\n" % item)


                            with open("Prod_envelope_Aero_Min"+i+".txt", 'w') as f:
                                for item in productMin:
                                    f.write("%s\n" % item)
                    #print(KI.keys())
                    #add each key to model
                    #KI_list = []
                    #for react_id in KI.keys():
                    #    KI_list.append(react_id)
                    #    r = universal.reactions.get_by_id(react_id)
                    #    model.add_reaction(r)
                    

                    
                    
                    # Production envelope model variant with only KIs
                    
                        #prod_env_KI = cobra.flux_analysis.phenotype_phase_plane.production_envelope(KI_model, biomass.id, objective='EX_'+target[0]+'_e', carbon_sources=c_source.id, points=10)
                        p_KI['only_KIs']=(growth, productMax, productMin)
                        prod_env_KI_KO.append(p_KI)
                    except Infeasible as inf:
                        print('''FVA is infeasible:''', inf)
                        p_KI['only_KIs']=('Infeasible')
                        prod_env_KI_KO.append(p_KI)


                    
                    #production envelope on different KOs strategies

                    print(values[i][0])
                    for n in range(len(values[i][0])):
                        knock_model = values[i][0][n]
                        biomass = get_biomass_equation(KI_model)
                        print('----'+str(n)+'-----')
                        p_KO = {} 
                        try:
                            growth=[]
                            productMax=[]
                            productMin=[]
                            #prod_env_KO = cobra.flux_analysis.phenotype_phase_plane.production_envelope(knock_model.model, biomass.id, objective='EX_'+target[0]+'_e', carbon_sources=c_source.id, points=10)
                            KOs = [r.id for r in knock_model.solution.KOs]
                            knock_model.model.objective=biomass 
                            #model.reactions.get_by_id(source[0]).lower_bound=-1 
                            maxG=knock_model.model.optimize().objective_value
                            knock_model.model.objective='EX_'+target[0]+'_e'
                            print(knock_model.model.objective.expression, target[0])
                            for o in range(Niter+1):
                                g=maxG/Niter*o
                                print(g)
                                biomass.lower_bound=g
                                maxpk = knock_model.model.optimize().objective_value
                                minpk = knock_model.model.optimize(objective_sense='minimize').objective_value
                                productMax.append(maxpk)
                                productMin.append(minpk)
                                growth.append(g)
                                
                                
                                with open("Prod_envelope_Aero"+i+"KO"+str(n+1)+".txt", 'w') as f:
                                    for item in growth:
                                        f.write("%s\n" % item)
                                
                        
                                with open("Prod_envelope_Aero_Max"+i+"_KO"+str(n+1)+".txt", 'w') as f:
                                    for item in productMax:
                                        f.write("%s\n" % item)
                                                

                                with open("Prod_envelope_Aero_Min"+i+"_KO"+str(n+1)+".txt", 'w') as f:
                                    for item in productMin:
                                        f.write("%s\n" % item)



                            p_KO[str(n+1)+'KO']=({'KOs': KOs, 'Objective': knock_model.prob.objective.value(), 'KO_model':knock_model}) #info_dictionary['KOs'] to get directly the KOs
                            prod_env_KI_KO.append(p_KO)
                            print('FVA worked')
                        except Infeasible as inf:
                            print('''FVA is infeasible:''', inf)
                            KOs = [r.id for r in knock_model.solution.KOs]
                            p_KO[str(n+1)+'KO']=({'KOs': KOs, 'Objective':knock_model.prob.objective.value(), 'KO_model':knock_model}, 'Infeasible')
                            prod_env_KI_KO.append(p_KO)
                    #remove_rlist(KI_list, model)
            production_envelope_per_model[i]=prod_env_KI_KO
    return production_envelope_per_model

def full_knock_out_analysis(input_file, output_consumption, output_con_and_prod, model, universal):
    """
    Find knockout strategies and evaluate them with production envelope

    combines the two functions run_optknock_analysis and 
    production_env_KO_eval which respectively run Optknock
    algorithm and evaluate the solutions by calculating the production
    envelope on the model variant without and with the reaction
    knock outs. 
    -----------------------------------------------------------------
    input_file--str input file in .csv format dictionary like
    output_consumption--dict as output from analysis_gf_sol function
    output_con_and_prod--dict retured by cons_prod_dict function
    model--cobra.Model reference model in BiGG namespace
    universal--cobra.Model universal model in BiGG namespace

    Return:
    evalKOs: dict with the pandas.DataFrame objects of the production 
        envelope of each model variants and the correspondant KO 
        mutants, as returned by production_env_KO_eval function
    """
    output_optknock = run_optknock_analysis(output_consumption, output_con_and_prod, model, universal)
    evalKOs = production_env_KO_eval(input_file, model, universal, output_optknock, output_con_and_prod)
    return evalKOs








	
