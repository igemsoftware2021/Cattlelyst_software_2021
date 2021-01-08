# Pipeline for the identification of metabolic engineering strategies and the suggestion of Biobrick parts

## Description
This pipeline is a Python package that allow users to fnd strain desing strategies in a genome scale metabolic model of interest. 

 The pipeline suggests engineering strategies to make a chassis organism able to **grow on atypical carbon sources** (e.g. pollutants such as methane) and to use them as the substrates for the **production of targets** (e.g. high-value compounds). 

### Abstract 
Biotechnology is particularly suited for targeting two societal needs: _the degradation of pollutants and the production of commercially relevant chemicals_. Microorganisms have great potential for the completion of both tasks. Their natural variability is harnessed by means of _metabolic engineering_ and _synthetic biology_ with the aim of creating strains with improved phenotypes for bio-degradation/production. Currently, numerous tools are available for the in silico design of microorganisms producing compounds relevant for the industry. Nonetheless, few tools exists for the identification of organism-specific strategies for strains to be used in _bioremediation_ Here we present a pipeline in Python that suggests engineering strategies allowing the **degradation of substrates such as pollutants and the simultaneous production of high-value compounds in genome scale metabolic models**. The engineering approaches are defined by studying models of metabolism using constraints-based methods. The resulting strategies account for both insertion and deletion of reactions. Moreover, the pipeline was built with the intent of facilitating the transition between the computational design and the validation of the strategies in vivo. For this purpose a wikibase of standardised biological parts was connected to the pipeline, which is thereby able to suggest the protein coding sequences that have to be introduced into the chassis for reproducing the designed pathway. Additionally, the Codon Harmonizer tool has been expanded and implemented in the pipeline with the purpose of generating codon-harmonised restriction sites-free Biobricks’s sequences. We validated the pipeline by using it for the identification of engineering designs in three test cases on the models of two organisms: Escherichia coli and *Pseudomonas putida*. The pipeline was also applied to an original engineering objective: the development of a strain of *Escherichia coli* utilizing methane as the sole carbon source and producing L-lactate. The results proved first of all that the pipeline can be used for the identification of candidate knock-ins and it is able to find unpredicted solutions. Ten engineering strategies were identified for L-lactate production in a methane-consuming E. coli strain. Secondly, the workflow for the evaluation of pathways’ thermodynamic feasibility was successfully tested. The pipeline we present provides an effective means for the rational design of chassis and their experimental testing in the future.

### Features
 
 The pipeline has three innovative aspects: 
 1. the prioritization of the assimilation of an uncommon substrate as the sole carbon source, which is particularly interesting for bioremediation purposes; 
 2. the identification of Biobrick parts corresponding to heterologous reactions and
 3. their codon harmonisation, which is a step usually kept separated from strain design tools. 
 
 The pipeline takes a user-filled template file indicating the specifics of the analysis as input.

### Flowchart
![Pipeline Workflow](/pipeline/workflow.png)

The pipeline's workflow is shown above.

- (i) Initialization involves loading the models, adding the exchange reactions and expanding the universal model; 
- (ii) Reaction search with the Gapfilling algorithm suggests reaction knock-ins. It is performed if the model does not already grow on the selected substrate and/or cannot produce the target compound; 
- (iii) Pathway definition and thermodynamic feasibility assessment evaluates the knock-in strategies by calculating the pathway’s Max-min driving force value [1]; 
- (iv) Comparison of the knock-in strategies is done by assigning scores with an additive function; 
- (v) Biobricks associated to the reactions knock-ins are identified querying the Biobrick wikibase [2];  
- (vi) The restriction site-free codon harmonised coding sequence of a Biobrick is obtained with an expanded version of Codon Harmonizer tool [3]; 
- (vii) Optionally knock-out strategies can be identified using Optknock [4]

## Usage
The pipeline can be used in Linux environment via command line call. Two arguments are taken: the path to the repository with the input data, which should correspond to ./pipeline/inputs, and the path to the .csv input file. Usage details can be displayed with ‘call_full_analysis.py --help’. 

_Example pipelie call from **Linux command line**_

<pre>
python3 call_full_analysis.py -r "./pipeline/inputs" -i "./pipeline/inputs/lac_201003_short.csv"
</pre>

The pipeline’s functions can be used from a Python IDE. Jupter Notebook is the one used for the test in pipeline/tests/ directory.

## Requirements
* Python version 3.6.10
* [COBRApy](https://github.com/opencobra/cobrapy/tree/devel/src/cobra) [5] version 0.18.1
* The solver used throughout the whole pipeline is GLPK [6], the one selected default by COBRApy, except for Optknock analysis, which requires **CPLEX solver** [7] version 12.10.0.0
* The scripts for querying the wikibase of Biobricks uses wikidataintegrator package 
* wikidateintegrator [8]
* [equilirator_pathway](https://gitlab.com/equilibrator/equilibrator-pathway) [9] v 0.3.1

## Status
The author of the pipeline is Delielena Poli

The pipeline is in development. The step v in the flowchart has to be finalized and implemented in the pipeline. 

## Similar tools
 
[cameo](http://cameo.bio/)

[OptStrain](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC525696/)

[OptCouple](https://www.sciencedirect.com/science/article/pii/S2214030118300373)

## References
[1] E. Noor, A. Bar-Even, A. Flamholz, E. Reznik, W. Liebermeister, and R. Milo, “Pathway Thermodynamics Highlights Kinetic Obstacles in Central Metabolism,” PLoS Comput. Biol., vol. 10, no. 2, p. e1003483, Feb. 2014, doi: 10.1371/journal.pcbi.1003483.

[2] I. Fundation, “iGEM foundation.” https://igem.org/Main_Page (accessed Oct. 24, 2020).

[3] N. J. Claassens et al., “Improving heterologous membrane protein production in Escherichia coli by combining transcriptional tuning and codon usage algorithms,” PLoS One, vol. 12, no. 9, Sep. 2017, doi: 10.1371/journal.pone.0184355

[4] A. P. Burgard, P. Pharkya, and C. D. Maranas, “OptKnock: A Bilevel Programming Framework for Identifying Gene Knockout Strategies for Microbial Strain Optimization,” 2003, doi: 10.1002/bit.10803.

[5] R. García-Granados, J. A. Lerma-Escalera, and J. R. Morones-Ramírez, “Metabolic Engineering and Synthetic Biology: Synergies, Future, and Challenges,” Front. Bioeng. Biotechnol., vol. 7, no. MAR, p. 36, Mar. 2019, doi: 10.3389/fbioe.2019.00036.

[6] GLPK, “GLPK.”Available: http://www.gnu.org/software/glpk.

[7] CPLEX.I.I.V12, “1:User’s Manual for CPLEX. International Business Machines Corporation,” vol. 46, no. 53, p. 157, 2009.

[8] A. Waa, “A protocol for adding knowledge to Wikidata, a case report,” bioRxiv, p. 2020.04.05.026336, Apr. 2020, doi: 10.1101/2020.04.05.026336.

[9] A. Flamholz, E. Noor, A. Bar-Even, and R. Milo, “eQuilibrator-the biochemical thermodynamics calculator,” Nucleic Acids Res., vol. 40, 2012, doi: 10.1093/nar/gkr874