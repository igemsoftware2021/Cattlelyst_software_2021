#!/usr/bin/env python

"""
Codon harmonizer
A program that attempts to recode a gene by considering the relative usage
of each codon in it's host and selecting a codon with the nearest relative
usage in a target organism.
In typical codon optimizers, each codon of a gene of interest is converted to
the "best" codon for a target organism. Yet, wild-type sequences don't always
use the "best" codon in their host organism. This code adjusts for this by
selecting the codon of a target organism that most closely approximates the 
codon's usage in a source organism.
Try `codonharmonizer.py --help` for usage details.
This code was intended as a fork of Bart Nijsse's `codon harmonizer`
(https://gitlab.com/wurssb/Codonharmonizer)
but ended up a whole-scale rewrite.
Where referenced in academic work, you may cite this repository and may also
consider referencing the manuscript discussing Nijsse's work
(doi.org/10.1371/journal.pone.0184355).
Author: Shyam Saladi (saladi@caltech.edu)
Date: May 2020
License: GPLv3
"""

import sys
import os.path
import argparse
import re
import textwrap
import itertools
import json
import warnings
import logging, sys

import pandas as pd

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(stream=sys.stderr, level=logging.ERROR)
re_atcg = re.compile("[ATUCG]", re.IGNORECASE)

#added by Delielena Poli in date 09/09/2020
RSs = {'EcoRI':'GAATTC', 'XbaI':'TCTAGA', 'PstI':'CTGCAG', 'SpeI':'ACTAGT'} 

def read_fasta_file(fasta_file):
    """Yield each dna sequence as pd.DataFrame of codons/aa"""
    for rec in Bio.SeqIO.parse(fasta_file, "fasta"):
        # Check for valid length
        if len(rec.seq) % 3 != 0:
            warnings.warn("Dropping sequences with incomplete codons. Check input")
            continue

        # Check for invalid DNA characters
        if len(re_atcg.sub('', str(rec.seq))) > 0:
            warnings.warn("Dropping sequences with non-ACTG characters. Check input")
            continue

        rec.seq = Seq(str(rec.seq).upper().replace('U', 'T'))
 
        yield rec

def format_fasta(name, seq, width=70):
    seq = textwrap.fill(seq, width)
    return ">{}\n{}".format(name, seq)

def write_reference_freq(fasta_file):
    """Processes fasta file into codon totals"""

    # Accumulate codon counts
    codon_count = {}
    for c in itertools.product('ACGT', repeat=3):
        c = "".join(c)
        codon_count[c] = 0

    for rec in read_fasta_file(fasta_file):
        for i in range(0, len(rec.seq), 3):
            c = str(rec.seq[i:i+3])
            codon_count[c] += 1

    codon_count = pd.Series(codon_count)
    codon_count.index.name = 'codon'
    codon_count.name = 'count'
    df = codon_count.reset_index()

    # Add amino acid to frame
    def translate_codon(x):
        rec = SeqRecord(seq=Seq(x))
        return str(rec.translate().seq)
    df['aa'] = df['codon'].apply(translate_codon)

    df = df[['codon', 'aa', 'count']]
    df.sort_values(['codon'], inplace=True)
    return df

def read_reference_freq(freq_file_path):
    """Relative frequencies and counts of all the codons from the given organism"""
    df_freq = pd.read_csv(freq_file_path)
    assert df_freq.columns.tolist() == ['codon', 'aa', 'count']
    
    df_freq['codon'] = df_freq['codon'].str.upper()
    df_freq['aa'] = df_freq['aa'].str.upper()
    df_freq['count'] = df_freq['count'].astype(int)

    # 'freq' is called the 'Relative Codon Adaptiveness (RCA) score' by some
    df_freq['freq'] = df_freq.groupby("aa")['count'].transform(lambda x: x/x.max())
    return df_freq

def harmonize_gene(seq, df_source_freq, df_target_freq, stats=False):
    """Harmonizes a DNA sequence given source and target codon frequencies"""

    # Convert to dataframe
    codons = [str(seq[i:i+3]) for i in range(0, len(seq), 3)]
    protseq = list(str(seq.translate()))
    df_seq = pd.DataFrame({'codon': codons, 'aa': protseq})

    # map source codons to target codons
    df_source_freq = df_source_freq.set_index('codon')
    df_target_freq = df_target_freq.set_index('codon')
    
    codon_trans = {}
    for codon, r in df_source_freq.iterrows():
        aa, codon_freq = r['aa'], r['freq']

        df_aa = df_target_freq[df_target_freq['aa'] == aa]
        codon_repl = (df_aa['freq'] - codon_freq).abs().idxmin() # check loop to add in here, before writing the dict 
        codon_trans[codon] = codon_repl
    logging.debug(codon_trans)
            
    df_seq['recode'] = df_seq['codon'].apply(lambda x: codon_trans[x]) 
    #logging.debug(df_seq)
    seq_harmonized = "".join(df_seq['recode'].tolist())

    if stats:
        df_seq['wt_source'] = df_source_freq.loc[df_seq['codon'], 'freq'].values
        df_seq['wt_target'] = df_target_freq.loc[df_seq['codon'], 'freq'].values
        df_seq['recode_target'] = df_target_freq.loc[df_seq['recode'], 'freq'].values

    return seq_harmonized, df_seq

def RS_check(seq, seq_harmonized, source_freq, target_freq): #function added by Delielena Poli in date 09/09/2020
    """
    Generate a restriction site-free harmonized sequence

    Given an harmonized sequence, it checks for the presence of 
    restriction sites and if any is found it substitute the first
    codon involved in the site with the codon with the next lower 
    impact on CHI difference. 
    -----------------------------------------------------------------
    Arguments:
    seq--str, native sequence
    seq_harmonized--str harmonized sequence
    source_freq--pd.DataFrame with the relative frequencies and counts
        of all the codons from the source organism
    target_freq--pd.DataFrame with the relative frequencies and counts
        of all the codons from the target organism

    Return:
    seq_harmo: harmonized sequences without restriction sites
    """
    # map source codons to target codons
    df_source_freq = source_freq.set_index('codon')
    df_target_freq = target_freq.set_index('codon')
    
    # define a dictionary with synonimous codons
    synonymous = {}
    for codon, r in df_source_freq.iterrows():
        aa, codon_freq = r['aa'], r['freq']

        df_aa = df_target_freq[df_target_freq['aa'] == aa]
        synonymous[codon] = df_aa
        
    result = substitute_codons_in_RS(seq_harmonized, df_source_freq, seq, synonymous)
    n = result[0]
    seq_harmo = result[1]
    
    return seq_harmo
    
def substitute_codons_in_RS(seq_harmonized, df_source_freq, seq, synonymous): #function added by Delielena Poli in date 09/09/2020
    """
    Cheks and eliminates restriction sites from an harmonized sequence

    Given an harmonized sequence, it checks for the presence of 
    restriction sites the 4 restriction sites (RSs) used in biobricks
    and if any is found it identifies the reading frame, hence the 
    codons involved in the RS. Consequently it sobstitute the 
    codon with a synonymous one with the next lower impact on the 
    codon harmonization index (CHI) difference. 
    -----------------------------------------------------------------
    Arguments:
    seq_harmonized--str harmonized sequence
    df_source_freq--pd.DataFrame with the relative frequencies and 
        counts of all the codons from the source organism with 
        'codon' column set as index
    seq--str, native sequence
    synonymous--dict with synonympus codons

    Returns:
    n: int, number of sustitution 
    RS_free_example or 
        seq_harmonized: str, harmonized sequences without restriction 
            sites. If modifications were needed RS_free_example is 
            returned otherwise seq_harmonized is returned.
    """
    RS_free_example = '' 
    n = 0 # number of substitutions
    for enzyme, r_site in RSs.items():
        if n == 0:
            RS_free_example = seq_harmonized
            logging.debug('N==0') 
            if r_site in RS_free_example:
                n+=1
                print('There is a RS {} for the enzyme {} at position {}'.format(r_site, enzyme, RS_free_example.index(r_site)))
                index = RS_free_example.index(r_site)
                logging.debug(index)
                remainder = index % 3
                #logging.debug('reminder = {}'.format(remainder))
                if remainder == 0:# the RS starts when a codon starts
                    triplet_h = RS_free_example[index:index+3]
                    logging.debug('divisible by 3')
                    triplet_native = seq[index:index+3]
                    #logging.debug(triplet_native)
                    for codon, r in df_source_freq.iterrows():
                        if codon == triplet_native:
                            #logging.debug(codon)
                            AA_source = r['aa']
                            freq_codon_source = r['freq']
                            #logging.debug(freq_codon_source)
                            df_aa = synonymous[triplet_native]
                            logging.debug(df_aa)
                            df_aa = df_aa.drop(triplet_h, axis=0)
                            logging.debug(df_aa)
                            codon_repl = (df_aa['freq'] - freq_codon_source).abs().idxmin()
                            logging.debug(codon_repl)
                            RS_free_example = RS_free_example[:index]+codon_repl+RS_free_example[index+3:]

                else:# the RS is internal to a codon
                    index_first_involved_codon = index - remainder
                    logging.debug('non div_by_3')
                    triplet_h = RS_free_example[index_first_involved_codon:index_first_involved_codon+3]
                    triplet_native = seq[index_first_involved_codon:index_first_involved_codon+3]
                    logging.debug(triplet_h)
                    logging.debug(triplet_native)
                    for codon, r in df_source_freq.iterrows():
                        if codon == triplet_native:
                            AA_source = r['aa']
                            freq_codon_source = r['freq']
                            df_aa = synonymous[triplet_native]
                            logging.debug(df_aa)
                            df_aa = df_aa.drop(index=triplet_h)
                            logging.debug(df_aa)
                            codon_repl = (df_aa['freq'] - freq_codon_source).abs().idxmin()
                            logging.debug(codon_repl)
                            RS_free_example = RS_free_example[:index_first_involved_codon]+codon_repl+RS_free_example[index_first_involved_codon+3:]
                            logging.debug(RS_free_example)

            else:
                    print('There is no restriction site for {} in the harmonized sequence. \nNumber of substitution = {}'.format(enzyme, n))
        elif n >= 1 and RS_free_example != '': # if there have already been substitutions uses the already modified sequence (RS_free_example) instead 
            logging.debug('N>= 1 ')
            logging.debug(RS_free_example)
            if r_site in RS_free_example:
                n+=1
                print('There is a RS {} for the enzyme {} at position {}'.format(r_site, enzyme, RS_free_example.index(r_site)))
                index = RS_free_example.index(r_site)
                logging.debug(index)
                remainder = index % 3
                logging.debug(remainder)
                if remainder == 0: # the RS starts when a codon starts
                    triplet_h = RS_free_example[index:index+3]
                    logging.debug(triplet_h)
                    triplet_native = seq[index:index+3]
                    logging.debug(triplet_native)
                    for codon, r in df_source_freq.iterrows():
                        if codon == triplet_native:
                            logging.debug(codon)
                            AA_source = r['aa']
                            freq_codon_source = r['freq']
                            logging.debug(freq_codon_source)
                            df_aa = synonymous[triplet_native]
                            logging.debug(df_aa)
                            df_aa = df_aa.drop(triplet_h, axis=0)
                            logging.debug(df_aa)
                            codon_repl = (df_aa['freq'] - freq_codon_source).abs().idxmin()
                            logging.debug(codon_repl)
                            RS_free_example = RS_free_example[:index]+codon_repl+RS_free_example[index+3:]

                else:# the RS is internal to a codon
                    index_first_involved_codon = index - remainder
                    logging.debug(index_first_involved_codon)
                    triplet_h = RS_free_example[index_first_involved_codon:index_first_involved_codon+3]
                    triplet_native = RS_free_example[index_first_involved_codon:index_first_involved_codon+3]
                    logging.debug(triplet_h)
                    logging.debug(triplet_native)
                    for codon, r in df_source_freq.iterrows():
                        if codon == triplet_native:
                            AA_source = r['aa']
                            freq_codon_source = r['freq']
                            df_aa = synonymous[triplet_native]
                            df_aa = df_aa.drop(triplet_h, axis=0)
                            logging.debug(df_aa)
                            codon_repl = (df_aa['freq'] - freq_codon_source).abs().idxmin()
                            logging.debug(codon_repl)
                            RS_free_example = RS_free_example[:index_first_involved_codon]+codon_repl+RS_free_example[index_first_involved_codon+3:]

            else:
                    print('There is no restriction site for {} in the harmonized sequence. \nNumber of substitution = {}'.format(enzyme, n))
            
    if RS_free_example != '':
        logging.debug(RS_free_example)
        return n, RS_free_example
    else:
        print('The sequence does not constain Restriction Sites')
        logging.debug(seq_harmonized)
        return n, seq_harmonized

def calculate_CHI(seq_target, seq_source, df_target_freq, df_source_freq):
    """Calculate Codon Harmonization Index (CHI)"""

    df_source_freq = df_source_freq.set_index('codon')
    df_target_freq = df_target_freq.set_index('codon')

    diffsum = 0
    for i in range(0, len(seq_source), 3):
        codon_target = str(seq_target[i:i+3])
        codon_source = str(seq_source[i:i+3])
        diffsum += abs(df_target_freq.loc[codon_target, 'freq'] -
                       df_source_freq.loc[codon_source, 'freq'])

    return diffsum/(len(seq_source)/3)

def main():
    parser = argparse.ArgumentParser(
        description='Harmonize a coding sequence for a target organism')
    parser.add_argument("fasta")

    parser.add_argument("--write_freqs", action='store_true',
        help="Count codons and write to stdout, e.g. of source or target coding sequences")

    arg_group = parser.add_argument_group('organism',
        "Coding sequences (fasta file) or codon counts (csv file generated by `--write_freqs`)")
    arg_group.add_argument("--source", required='--write_freqs' not in sys.argv)
    arg_group.add_argument("--target", required='--write_freqs' not in sys.argv)

    parser.add_argument("--stats", default=None,
        help="Write statistics for harmonized genes as a json-formatted file")

    args = parser.parse_args()

    if args.write_freqs:
        df = write_reference_freq(args.fasta)
        print(df.to_csv(index=False))
        return

    # Read source and target frequencies
    source_name, _ = os.path.splitext(args.source)
    source_freq = read_reference_freq(args.source)
    target_name, _ = os.path.splitext(args.target)
    target_freq = read_reference_freq(args.target)

    if args.stats:
        stats = []
    
    # harmonize all sequences provided
    for rec in read_fasta_file(args.fasta):
        header, seq = rec.id, rec.seq
        seq_harm, df_stats = harmonize_gene(seq, source_freq, target_freq, args.stats)
        # statistics
        initial_CHI = calculate_CHI(seq, seq, target_freq, source_freq)
        harm_CHI = calculate_CHI(seq_harm, seq, target_freq, source_freq)

        name = "{} (src:{}-tgt:{}) (CHI_0:{}-CHI:{})".format(
            header, source_name, target_name, initial_CHI, harm_CHI)
        print(format_fasta(name, seq_harm))
        
        if args.stats:
            with open(args.stats, 'a') as fh:
                json.dump({
                    'name': header,
                    'sequences': {
                        'initial': str(rec.seq),
                        'harmonized': seq_harm,
                    },
                    'CHI': {
                        'initial': initial_CHI,
                        'harmonized': harm_CHI,
                    },
                    'data': {c:df_stats[c].tolist() for c in df_stats},
                }, fh)
    
    return

if __name__ == "__main__":
    main()
