# Description of scripts 
These scripts constitute the pipeline's modules

The modules used by the pipeline consist in eleven Python files that are contained in the “scripts” directory. Three of those files are scripts that have been only partly modified:
- Optknock_robustknock.py
- call_Optknock_robustknock.py
- codon_harmonizer_RSs.py

The rest have been originally written for the pipeline. The script for querying the Biobrick wikibase still has to be added and will be done so when wikibase development is complete. 

The single script that connects all the modules together is call_full_analysis.py in the main repository of the project (outside scripts directory). That script is the one that can be calle via command line.