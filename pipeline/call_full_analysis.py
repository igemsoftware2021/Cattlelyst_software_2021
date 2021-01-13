#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script calling the pipeline from command line  

Find reactions KIs (and KOs) allowing growth on uncommon substrate
(and production of a target compound).
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"


import csv
import argparse
import sys
import cobra
from scripts.analysis import analysis_gf_sol, dict_prod_sol, cons_prod_dict
from scripts.import_models import get_universal_main, get_reference_model
from scripts.input_parser import main
from scripts.scoring_system import (scores_evaluations, score_BB_presence, 
                            plotting_KI_BB_scores)
from scripts.to_and_from_biobrick_wikidata import (output_for_biobrick_search, \
                                        run_harmonization_per_model)
from scripts.call_Optknock import (run_optknock_analysis, \
                                        add_Optknock_to_analysis)
from scripts.output_generation import (summary_output, metab_out_chain, 
                            write_data_provenance)

from query_bioparts import query_bioparts

def parse_options():
    """
    Define arguments for command line call of the whole analysis

    Uses argparse package to define the two arguments for the analysys:
    -r the respoitory where the models (reference and universal) are
        stored
    -i the name (and path) to the inpur file in .csv format.
    """
    usage="\ncall_full_analysis.py -r <path to repository> -i <input csv file>"
    description="Find reactions KIs (and KOs) allowing growth on uncommon substrate and production of a target compound."
    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument("-r", "--repository", dest="data_repo",
                    help="string of the path to the repository with the models", required=True)
    parser.add_argument("-i", "--input_file", dest="input_file",
                    help="string of the input file name with .csv extension specified", required=True)
    args = parser.parse_args()

    return args

def main_whole_analysis(argv):
    """
    whole analysis workflow.

    Calls all the functions for the whole analysis:
    1. load the models
    2. read the import file and correct the model accordingly
    3. reaction finding with gapfilling to allow growth on desired 
        compound
    4. reaction finding with gapfilling to allow production of the 
        target 
    5. evaluate the solutions 
    6. create input for querying biobrick wikidata
    7. harmonize the sequences of the matching biobricks
    8. run optknock analysis if it is selected by the user

    It returns output files, saved in the path "./outputs/". The two
    final files are:
    - summary_output.csv
    - detailed_output.csv
    - data_provenance.txt
    Other intermediate files such as consumption.csv can also be found
    at the same location. 
    ------------------------------------------------------------------
    Arguments:
    argv--arguments from argparse functions and interaction with sys
        module
    """
    inputs=parse_options()
    model = get_reference_model(inputs.data_repo, inputs.input_file)
    universal = get_universal_main(inputs.data_repo, inputs.input_file)

    main(inputs.input_file, universal, model)

    consumption = analysis_gf_sol(inputs.input_file, model, universal)

    model.remove_reactions(['LCADi']) # TODO: delete this line in final version

    production = dict_prod_sol(inputs.input_file, consumption, model, universal)

    # Get overall solutions from GapFilling analysis, hence the 
    # different strategies
    final = cons_prod_dict(inputs.input_file, model, universal, consumption, production)
    # Score the metabolic engineering strategies
    scores_output = scores_evaluations(inputs.input_file, consumption, final, optknock_analysis=False)
    # Formatting information for wikidata query   
    to_wikidata = output_for_biobrick_search(universal, scores_output, final)
    # Retrieve from_wiki using query_bioparts function
    from_wiki = query_bioparts(to_wikidata)

    # Score Biobrick presence
    ##scores_BB = score_BB_presence(inputs.input_file, to_wikidata, from_wiki, scores_output)
    # Plot scores for knock-ins engineering strategy 
    # vs scores for presence of matching biobrick 
    ##plotting_KI_BB_scores(scores_output, scores_BB)
    # Write the scores to the summary output csv file
    ##summary_output(scores_output, scores_BB)
    # Harmonize the non native sequences per model variant (hence per engineering strategy)
    ##models_harmonized_seq = run_harmonization_per_model(from_wiki, inputs.input_file)     
    ##print(models_harmonized_seq)
    # Generate detailed output.   
    ##metab_out_chain(inputs.input_file, consumption, final, scores_output, to_wikidata, models_harmonized_seq)
    

    # Read the input file to decide whether to perform Optknock or not
    consider_Optknock = add_Optknock_to_analysis(inputs.input_file)

    write_data_provenance(inputs.data_repo, inputs.input_file)

    if consider_Optknock:
        
        print('Following OptKnock analysis...')
        result_optknock = full_knock_out_analysis(consumption, final, model, universal)
        # TODO: plot the results
        summary_output(scores_output, scores_BB, out_opt=result_optknock)
        

if __name__ == "__main__": 
    main_whole_analysis(sys.argv[1:])

        


