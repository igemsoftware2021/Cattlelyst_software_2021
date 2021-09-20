#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module to generate output files.

This modules writes the csv files composing the final
output of the pipeline.
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"

import logging
import sys

import csv
import cobra
import pandas as pd
import platform
import pkg_resources
from datetime import datetime

from pipeline.scripts.analysis import get_c_source

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)

def summary_output(scores_output, scores_bb, out_opt=None):
    """
    Create summary_output.csv file

    It overwrite the file generate with generate_scores function
    and adds the rank to it. 
    -------------------------------------------------------
    Arguments:
    scores_output--dict returned by scores_evaluation function
    scores_bb--dict retured by score_BB_presence function
    opt--bool if True the summary will contain information 
        from Optknock analysis

    Optional
    out_opt--dict returned by run_optknock_analysis function
        default None
    """
    #transform the dictionary from score evaluation function to 
    # pandas.DataFrame and add it to csv file.
    new_df = pd.DataFrame.from_dict(scores_output)
    rank = new_df.transpose()
    rank.to_csv("pipeline/outputs/partof_summary_output.csv", index = True, mode='a', header=False)
    
    #read the full csv file as a DataFrame
    df = pd.read_csv("pipeline/outputs/partof_summary_output.csv")
        
    #add the percentage of biobricks to df  
    for ind, df_ind in enumerate(df.loc[df['tot_metabolic'].notnull()].index.values):
        for tup in scores_bb:
            if tup[0] == df.iloc[df_ind,0]:
                df.at[df_ind,'percentage_Biobrick'] = tup[1]

    if out_opt!=None:
        #add info from optknock to df
        for i, values in out_opt.items():
            for inx, df_ind in enumerate(df.loc[df['model_name'].notnull()].index.values):
                if i == df.iloc[df_ind,0]:
                    for index, dic in enumerate(values):
                        for n in range(1,4):
                            if index==n:
                                obj = dic[str(n)+'KO'][0]['Objective']
                                all_ko = str(n)+'KO'
                                r_kos = dic[str(n)+'KO'][0]['KOs']
                                print(all_ko, r_kos, obj)
                                df['allowed_KOs'][df_ind+n-1] = all_ko
                                df['reactions_KOs'][df_ind+n-1] = r_kos
                                df['objective_value'][df_ind+n-1] = obj

    df.to_csv("pipeline/outputs/summary_output.csv", index = False, mode='a', header=True)

def metabolic_detailed_output(input_file, consumption, final):
    """
    Write information on reaction addition to the detailed output file

    Read the dictionaries from consumption and and the overal reaction 
    addition analysis and it will generate a file with the reaction ID 
    and rates.

    It fills the following columns of detailed_output.csv file:
    - model_name,	
    - uptake_rate,	
    - production_rate
    - MDF 
    - pathway_length	
    - reactions_added	
    - fluxes
    -------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    consumption--dict retured by analysis_gf_sol function
    final--dict retured by cons_prod_dict function
    """
    c_source = get_c_source(input_file)
    with open('pipeline/outputs/metabolic_detailed_output.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'model_name', 'uptake_rate', 'production_rate',
            'MDF', 'pathway_length', 'reactions_added', 'fluxes', 
            'EC_number', 'KEGG_ID', 'Biobrick_ID', 
            'Biobrick_sequence', 'Native_organism',
            'Harmonized_seq', 'CHI_diffeence'
            ])
        writer.writerow([
            "#identifier of the model: the first number indicate the variant consumption's reaction search, this is followed by the identifier of the target, finally, the last number indicates the variant for production's reaction search", '#flux through the exchange reactions of the carbon source. Unit = mmol/gDW/h. Negative values indicate uptake', '#flux throught the exchange reaction of the targetUnit = mmol/gDW/h. Positive values indicate production.',
            '#MDF value. Unit = kJ/mol. The higer the MDF value, the higher the thermodynamic driving force', '#number of reactions steps from substrate to product. The lower the better', '#identifiers of the BiGG reactions found with gapfilling. They represent the reactions knock ins (KIs) suggested to reach growth on the selected c-source and production of the target', '#fluxes through the reactions KIs. Unit = mmol/gDW/h ', 
            '#EC number of the reactions KIs', '#KEGG reactions IDs of the KIs', '#identifier of the Biobrick that carries gene(s) encoding for enzyme(s) whose functions matches the reactions KIs', 
            '#coding sequence in fasta format of the gene(s) included in the biobrick part with the matching function/reaction', '#species name of the organis from which the gene has been cloned',
            "#codon harmonized gene sequence. The sequence provided in 'Biobrick sequence' column is codon ahrmonized and it is free from restriction sites used in biobrick assembly (EcoRI, SpeI, XbaI, PstI)", '#native codon harminizaition index (CHI) - codon harmonized and restriction site free CHI. The closer to 0 the difference is, the best the codon harmonization with the native frequences is'
            ])
        for n in range(1, len(consumption)+1):
            #print('---------'+str(n)+'---------')
            up_rate = [consumption[n][1][3]['EX_'+c_source+'_e']]
            for key, values in final['production_'+str(n)].items():
                p_rate = ''
                x = values[3].values()
                for k in x:
                    p_rate += str(k)
                r_added = []
                fluxes = []
                dict_r_added = values[4]
                for r, f in dict_r_added.items():
                    r_added.append(r)
                    fluxes.append(f)

                thermo = values[5]
                for m in range(len(r_added)):
                    if (m+1) <= len(up_rate):
                        writer.writerow([
                        (str(n)+'_'+key), consumption[n][1][3]['EX_'+c_source+'_e'], float(p_rate), 
                        thermo['mdf'], thermo['pathway_length'], r_added[m], fluxes[m]
                        ])
                    else:
                        writer.writerow(['', '', '', '', '', r_added[m], fluxes[m]])

        writer.writerow(['---------------------------------------------'])    
     

def reaction_details(scores_output, final, to_wiki, out_harm):
    """
    Add information of reactions and associated biobricks and sequences

    Reads the intermediate file metabolic_detailed_output.csv as 
    pandas.DataFrame and per each reaction identifiers reads the
    corresponding EC number and KEGG IDs from to_wiki dictionary. 
    Then, it reads the output of the sequence harminization and adds 
    information on the biobricks matching each EC number.

    The DataFrame with the additional information is saved as the final 
    detailed_output.csv file 
    -------------------------------------------------------
    Arguments:
    scores_output--dict returned by scores_evaluation function
    final--dict retured by cons_prod_dict function
    to_wiki--list Nested list with information of EC and KEGG 
        per reaction.returned by output_for_biobrick_search function
        from to_and_from_biobrick_wikidata module. 
    out_harm--nested dictionary resulting from run_harmonization_per_model
        function. The model variants as keys, and the EC numbers of the 
        heterologous reactions as keys of the nested dictionary, while 
        the values are tuples.
    """
    df = pd.read_csv("pipeline/outputs/metabolic_detailed_output.csv")
    df.loc[:, 'EC_number'].astype("string")
    for n in range(1,len(scores_output)+1):
        logging.debug('n = {}'.format(n))
        info_for_model=[]
        solution = scores_output[n]
        vals = to_wiki[n-1]
        added_reactions = final['production_'+str(solution[0][0])][str(solution[0][2:])][4].keys()
        logging.debug(added_reactions)
        for inx, df_ind in enumerate(df.loc[df['model_name'].notnull()].index.values):
            start = df_ind
            if solution[0] == df.iloc[df_ind,0]:
                for ind, i in enumerate(added_reactions):
                    if ind == 0:
                            logging.debug(ind)
                            info = vals[ind]['reaction']
                            logging.debug(info)
                            if info != None:
                                logging.debug(solution[0])
                                logging.debug(df.iloc[df_ind,0])
                                logging.debug(df.loc[df_ind,'model_name'])
                                df['EC_number'][df_ind] = str(info['EC'])
                                df['KEGG_ID'][df_ind] = str(info['KEGG'])

                            else:
                                df['EC_number'][df_ind] = 'N/A'
                                df['KEGG_ID'][df_ind] = 'N/A'

                    else:
                        for k in range(1,len(added_reactions)):
                            info = vals[k]['reaction']
                            if info != None:
                                logging.debug(solution[0])
                                df['EC_number'][df_ind+k] = str(info['EC'])
                                df['KEGG_ID'][df_ind+k] = str(info['KEGG'])

                            else:
                                df['EC_number'][df_ind+k] = 'N/A'
                                df['KEGG_ID'][df_ind+k] = 'N/A'
    for _ in range(len(scores_output)):
        for key, values in out_harm.items():
            for i, v in values.items():
                for inx, df_ind in enumerate(df.loc[df['EC_number'].notnull()].index.values):
                    if i == df.loc[df_ind, 'EC_number']:
                        if v!=None:
                            df['Biobrick_ID'][df_ind] = str(v[0])
                            df['Biobrick_sequence'][df_ind] = str(v[1])
                            df['Native_organism'][df_ind] = str(v[2])
                            df['Harmonized_seq'][df_ind] = str(v[3])
                            df['CHI_diffeence'][df_ind] = str(v[4])

                        else:
                            df['Biobrick_ID'][df_ind] = 'N/A'

    df.to_csv("pipeline/outputs/detailed_output.csv", index = False, mode='a', header=True)

def metab_out_chain(input_file, consumption, final, scores_output, to_wiki, out_harm):
    """
    Fill detailed_output.csv

    Call of metabolic_detailed_output and reaction_details to generate
    Fill the csv file with reaction detailes (EC and KEGG ID) per
    model.
    -------------------------------------------------------
    Arguments:
    input_file--str input file in .csv format dictionary like
    consumption--dict retured by analysis_gf_sol function
    final--dict retured by cons_prod_dict function
    scores_output--dict returned by scores_evaluation function
    to_wiki--list nested list with information of EC and KEGG per
        reaction.
    out_harm--nested dictionary resulting from run_harmonization_per_model
        function. The model variants as keys, and the EC numbers of the 
            heterologous reactions as keys of the nested dictionary, while 
            the values are tuples.
    """
    metabolic_detailed_output(input_file, consumption, final)
    reaction_details(scores_output, final, to_wiki, out_harm)

def write_data_provenance(data_repo, input_file):
    """
    Generate .txt file for data provenance of inputs and outputs 

    This function is called every time the pipeline is run, and
    it creates a .txt file with information on time and date, inputs,
    operating system, python and cobrapy version and outputs.
    """
    with open('pipeline/outputs/data_provenance.txt', 'w') as file:
        #write time and date
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        file.write('Pipeline run at the following date and time:\n')
        file.write(dt_string)
        file.write('\n')
        #info on input
        file.write('\nThe input was:\n')
        file.write(input_file)
        file.write('\nPath to file:\n')
        path = data_repo+input_file
        file.write(path)
        file.write('\n')
        #write platfrom info
        file.write('\nSystem network name:\n')
        file.write(platform.node())
        file.write('\nPlatform information:\n')
        file.write(platform.platform())
        file.write('\n')
        #write python version
        file.write("\nPython version:\n")
        file.write(sys.version)
        #write cobrapy version
        cobra_version = pkg_resources.require("cobra")[0].version
        file.write("\nCobrapy version:\n")
        file.write(cobra_version)
        file.write('\n')
        #output location
        file.write("\nAside from this data_provenance.txt file,\n")
        file.write("The following output files can be found \n")
        file.write("in the current repository (pipeline/ouptuts/<filename>):")
        file.write("""\n
            - summary_output.csv
            - detailed_output.csv
            - consumption.csv
            - metabolic_detailed_output.csv
            - partof_summary_output.csv\n""") 