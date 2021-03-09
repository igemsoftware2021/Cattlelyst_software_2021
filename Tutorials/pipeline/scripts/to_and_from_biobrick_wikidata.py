#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for interacition with Biobrick wikidata page and output evaluation

Contains functions to query the biobrick wikidata and parse the output.
The output is used to:
1) get the percentage of biobrick matching the heterologous reactions 
    of a pathway 
2) generate codon harmonized and restriction sites free sequences of the
    gene(s) in the matching biobricks.
"""
__author__ = "Delielena Poli"
__email__ = "delielena.poli@wur.nl"
__status__ = "Development"

import cobra
import csv
import warnings
import gzip
import os
import wget
import logging, sys

from pipeline.scripts.import_models import get_expression_host
from pipeline.scripts.codon_harmonizer_RSs import (read_reference_freq,
                            read_fasta_file, harmonize_gene,
                            RS_check, calculate_CHI, format_fasta)

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)

# To biobrick wikidata 

def output_for_biobrick_search(universal, scores, output_con_and_prod):
    """
    Create output for the query of the biobrick wikidata

    Read the results of the metabolic optimization and create output 
    for the query of the biobrick wikidata. The information needed 
    for the query is the EC number of the heterogous reactions found
    with gapfilling optimization.
    -----------------------------------------------------------------
    Arguments:
    universal--cobra.Model universal model in BiGG namespace
    scores--dict with the ranked model variants and associated scores. 
    output_con_and_prod--dict with comprehensive information on 
        reactions additions and fluxes for each model veraint 
        growing on the indicated substrate and producing the target.

    Return:  
    for_biobrick_search: list of lists with information on EC number
        and KEGG ID of each heterologous reaction
    """
    for_biobrick_search=[]
    for n in range(1,len(scores)+1):
        info_for_model=[]
        solution = scores[n]
        logging.debug(solution[0])
        added_reactions = output_con_and_prod['production_'+str(solution[0][0])][str(solution[0][2:])][4].keys()
        logging.debug(added_reactions)        
        for i in added_reactions:
            info_for_BB = {}
            # retrieve reactions from universal
            logging.debug(i)
            reaction = universal.reactions.get_by_id(i)
            if 'KEGG Reaction' in reaction.notes.keys() and 'EC Number' in reaction.notes.keys():
                keggid = reaction.notes['KEGG Reaction']
                ec = reaction.notes['EC Number']
                info_for_BB['reaction']={'EC':ec, 'KEGG':keggid}
            elif 'KEGG Reaction' not in reaction.notes.keys() and 'EC Number' in reaction.notes.keys():
                ec = reaction.notes['EC Number']
                info_for_BB['reaction']={'EC':ec, 'KEGG':None}
            elif 'KEGG Reaction' in reaction.notes.keys() and 'EC Number' not in reaction.notes.keys():
                keggid = reaction.notes['KEGG Reaction']
                info_for_BB['reaction']={'EC':None, 'KEGG':keggid}
            else:
                info_for_BB['reaction']= None
                logging.debug('NONE!')
                logging.debug(i)
            info_for_model.append(info_for_BB)
        for_biobrick_search.append(info_for_model)
    return for_biobrick_search

# From Biobrick Wikidata if doesn't work look 2020-09-15 connect jupy

def get_entries_for_species(spec):
    """
    Get the entries for each species from NCBI genome report
    
    It creates a list of dictionaries with the info on the entries 
    matching the species name in the prokaryote.txt file from NCBI.
    The file can be downloaded from this link: 
    ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
    -----------------------------------------------------------------
    Argument:
        
    spec--string, obtained from querying the wikidata

    Return:
    species_list: list of dictionaries. It is created for each heterologous 
            reaction and it contains the following keys:
            - species: name of the organisms in the entry
            - date: last modification date
            - status: type of genome assembly (Complete Genome, 
            scaffold, Contig)
            - assembly: FTP path to the ncbi assembly page
    
    """
    species_list = []

    with open('prokaryotes.txt', newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel', delimiter='\t')
        for row in reader:
            if spec != None and spec in row['#Organism/Name']:
                raw_list ={}
                raw_list['species'] = row['#Organism/Name']
                raw_list['date'] = row['Modify Date']
                raw_list['status'] = row['Status']
                raw_list['assembly'] = row['FTP Path']
                species_list.append(raw_list)
            else:
                continue
    #TODO: delete the following part for eukaryotes
    #with open('eukaryotes.txt', newline='', encoding='utf-8-sig') as csvfile:
    #    reader = csv.DictReader(csvfile, dialect='excel', delimiter='\t')
    #    for row in reader:
    #        if spec != None and spec in row['#Organism/Name']:
    #            raw_list ={}
    #            raw_list['species'] = row['#Organism/Name']
    #            raw_list['date'] = row['Modify Date']
    #            raw_list['status'] = row['Status']
    #            access = row['Assembly Accession']
    #            ftp_path = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/'+access[:3]+'/'+access[4:7]+'/'+access[7:10]+'/'+access[10:13]+'/'+access
    #            raw_list['assembly'] = ftp_path
    #            species_list.append(raw_list)
    #        else:
    #            continue

    return species_list

def get_link_to_assembly(species_list):
    """
    Get the path for assembly download from NCBI

    Select the most up-to-date Complete Genome of the each organisms
    and returns the FTP path for that one.
    -----------------------------------------------------------------
    Argument:

    species_list--list of dictionary obtained from 
        get_entries_for_species function
    
    Return:
    link_assembly: str, FTP path to the genome assembly page of the 
        selected model:
        e.g ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/448/615/GCA_010448615.1_ASM1044861v1
    """
    
    genomes = []
    link_assembly = ''
    for i in species_list:
        if i['status'] == 'Complete Genome' and i['assembly'] not in ('', '-'):
            genomes.append(i)

    date_list = []
    if len(genomes) != 0:
        for i in genomes:
            if len(genomes) == 1:
                link_assembly += i['assembly']
            elif len(genomes) >1:
                date_list.append(i['date'])
    else: 
        # consider eliminating the options below: if there is no 
        # genomic information, there isn't any cds_from_genomic file
        warnings.warn('There is no complete genome for this species')

    if len(date_list) != 0:
        latest = max(date_list)
        for i in species_list:
            if i['date'] == latest:
                link_assembly += i['assembly']
                break # to pick the first match in case more than one entry have the same date
        
    
    return link_assembly

def write_correct_link_to_cds(raw_link):
    """
    Get full link to specific NCBI CDS from gemonic file

    Gives the full link to the CDS from genimc data. This is the type 
    of input needed for calculating the codon frequencies used by codon 
    harmionizer. 

    This function is usefull mostly if the raw_link correspont to a 
    complete genome, since scaffold and contigs entries do not 
    necessarily have the file with extension _cds_from_genomic.fna.gz
    --------------------------------------------------
    Argument:
    raw_link--str, FTP path obtained from get_link_to_assembly function.

    Return:
    correct_link: it adds the extension _cds_from_genomic.fna.gz to the
                    raw_link. 
    e.g ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/448/615/GCA_010448615.1_ASM1044861v1/GCA_010448615.1_ASM1044861v1_cds_from_genomic.fna.gz
    """
    
    ID = raw_link.split('/')[-1]
    to_cds = ID + '_cds_from_genomic.fna.gz'
    
    correct_link = raw_link+'/'+to_cds
    return correct_link 

def get_native_info(spec):
    """
    Get correct link for each species

    Combins the three functions get_entries_for_species, 
    get_link_to_assembly and write_correct_link_to_cds
    -----------------------------------------------------------------
    Argument:
    spec--str, obtained from querying the wikidata

    Return:
    full_link: link to file with _cds_from_genomic.fna.gz extension
    """
    entries_list = get_entries_for_species(spec)
    #logging.debug(entries_list)
    #if len(entries_list) != 0 
    link = get_link_to_assembly(entries_list)
    if link != '' and link != '-' and link != '/':
        full_link = write_correct_link_to_cds(link)
        return full_link
    else:
        return None

def read_gene_seq(subdict):
    """
    Read gene sequence (str) and write it into a txt file 

    Given the output of the biobrick search by querying the wikidata,
    this function reads the string sequence of the gene encoding the
    enzyme with the desired function. It transfroms this sequence into
    a .txt file, for compatibility with the existing codon_harmonizer 
    workflow. 

    The file name has .fasta extension and is called sequence.fasta. 
    The header in the fasta format is >gene_seq
    -----------------------------------------------------------------
    Argument:
    subdict--dict, subdict from the output returned from wikidata query
        function. 
    """
    string = subdict['BB_Seq']
    with open("sequence.fasta", "w") as text_file:
            text_file.write(">gene_seq\n")
            text_file.write(string)

def get_gz_file_expression_host(input_file):
    """
    Read path to CDS from genomic of expression host

    Parse the input file to get the path to the file with the
    CDS of the expression host 
    -----------------------------------------------------------------
    Argument: 
    input_file--str, path to the fast .gz file with the CDS of the expression
        host
    Return: 
    file_name: str, name of the .gz fasta file as introduced by the user 
    """
    with open(input_file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel')
        for row in reader:
            if row['path_to_fna_file_expression_host'] != '':
                file_name = row['path_to_fna_file_expression_host']
    return file_name

def gunzip(source_filepath, dest_filepath, block_size=65536):
    """
    Unzip fasta.gz file from expression host

    Unzip the file to get access to CDS from genomic of the expression
    host. Unpacks and writes the file in the desired location
    ------------------------------------------------------------------
    Arguments:
    source_filepath--str, file with fna.gz extension 
    dest_filepath--str, path and file name of the .fna CDS file
        without .gz extensnion    
    """
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)

def run_harmonization_per_model(output_wikidata, input_file):
    """
    Read wikidata output and add harmonized sequences
    
    For each model and each EC number (which are associated to 
    the reaction) it obtains the coding sequences of native (with wget 
    package) and expression host and it derives the codon frequences.
    
    The codon harmonization is done on the fasta sequence of the 
    gene(s) encoded by the matching biobrick. The sequence is also 
    provided in the wikidata output.

    It returns a dictionary with the codon-harmonized and 
    restriction-free sequences for each model.
    -----------------------------------------------------------------
    Arguments: 
    output_wikidata--dict type (nested dictionary) returned from 
        wikidata query function. 
    input_file--str input file in .csv format dictionary like

    Return:
    harmonized_seq_per_model: nested dictionary with the model variants
        as keys, and the EC numbers of the heterologous reactions as 
        keys of the nested dictionary, while the values are tuples.            
    
    """
    n=0
    list_for_models = {}
    list_links = []
    harmonized_seq_per_model ={}
    exp_host = get_gz_file_expression_host(input_file)
    spec_exp_host = get_expression_host(input_file)
    #model_exp_host = get_ID_reference_model(input_file)
    for i in output_wikidata:
        n+=1
        list_for_models['model'] = n
        dict_links = {}
        harmonized_seq_per_reaction = {}
        for key, subdict in i.items():
            native = subdict['species']
            logging.debug(native)
            download_link = get_native_info(native)
            logging.debug(download_link)
            dict_links[native] = download_link
            if download_link != None:
                # do condon harmonization 
                # 1st get the CDS
                ## Native
                print('Beginning file download with wget module')
                reformatted = native.replace(' ','_')
                file_name = reformatted + '.cds.fna.gz'
                wget.download(download_link, file_name)
                new_file_name = file_name[0:-3]
                gunzip(file_name, new_file_name, block_size=65536)
                csv_file_name = reformatted+'.freq.csv'
                message = "codonharmonizer {} --write_freqs > {}".format(new_file_name, csv_file_name)
                logging.debug(message)
                os.system(message)
                ## Expression host
                gunzip(exp_host, exp_host[:-3], block_size=65536)
                fast_file_cds_exp_host = exp_host[:-3]
                os.system("codonharmonizer {} --write_freqs > host.freq.csv".format(fast_file_cds_exp_host))
                #get_frequences
                source_freq = read_reference_freq('native.freq.csv')
                target_freq = read_reference_freq('host.freq.csv')

                # read gene sequence
                read_gene_seq(subdict)
                # run codon harmonization
                

                for rec in read_fasta_file('sequence.fasta'):
                    header, seq = rec.id, rec.seq
                    seq_harm, df_stats = harmonize_gene(seq, source_freq, target_freq)
                    #check RS presencs
                    checked_harm_seq = RS_check(seq, seq_harm, source_freq, target_freq)
                    # add result of CH to a dict per model (and per reaction?)
                    initial_CHI = calculate_CHI(seq, seq, target_freq, source_freq)
                    harm_CHI = calculate_CHI(seq_harm, seq, target_freq, source_freq)
                    harm_CHI_RS_free = calculate_CHI(checked_harm_seq, seq, target_freq, source_freq)

                    name = "{} (src:{}-tgt:{}) (CHI_0={}-CHI={})".format(
                        header, native, spec_exp_host, initial_CHI, harm_CHI)
                    print(format_fasta(name, seq_harm))
                    RS_free = "{} (src:{}-tgt:{}) (CHI_0={}-CHI={})".format(
                        header, native, spec_exp_host, initial_CHI, harm_CHI_RS_free)
                    print(format_fasta(RS_free ,checked_harm_seq))
                    harmonized_seq_per_reaction[key] = (subdict['BB_name'], 
                                                        subdict['BB_Seq'], 
                                                        native, 
                                                        checked_harm_seq, 
                                                        'CHI_native={}-CHI_harmonized={}'.format(initial_CHI, 
                                                                                       harm_CHI_RS_free)
                                                       )
            else:
                harmonized_seq_per_reaction[key] = subdict['BB_Seq']
        harmonized_seq_per_model['Model {}'.format(n)] = harmonized_seq_per_reaction


        list_links.append(list_for_models)
        list_links.append(dict_links)
    return harmonized_seq_per_model