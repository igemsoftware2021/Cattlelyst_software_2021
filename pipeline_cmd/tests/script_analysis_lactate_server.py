#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for pipeline analysis in Linus server 

Find reactions KIs allowing for the test case on E. coli 
growing on methane and producing L-lactate.
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"

import csv

import cobra

from scripts.analysis import analysis_gf_sol, dict_prod_sol, cons_prod_dict
from scripts.import_models import get_universal_main, get_reference_model
from scripts.input_parser import main
from scripts.scoring_system import scores_evaluations
from scripts.to_and_from_biobrick_wikidata import output_for_biobrick_search


if __name__ == "__main__":

    model = get_reference_model('pipeline/inputs/', 'pipeline/inputs/lac_200925.csv')

    universal = get_universal_main('pipeline/inputs/', 'pipeline/inputs/lac_200925.csv')

    main('pipeline/inputs/lac_200925.csv', universal, model)

    consumption = analysis_gf_sol('pipeline/inputs/lac_200925.csv', model, universal, 'output_lac_4_iter.csv')

    print(consumption)

    model.remove_reactions(['LCADi'])

    production = dict_prod_sol('pipeline/inputs/lac_200925.csv', consumption, model, universal)

    print(production)

    final = cons_prod_dict('pipeline/inputs/lac_200925.csv', model, universal, consumption, production)

    scores_output = scores_evaluations('pipeline/inputs/lac_200925.csv', consumption, final)

    print(scores_output)

    to_wikidata = output_for_biobrick_search(universal, scores_output, final)

    


