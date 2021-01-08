# output directory

## Content description
The “outputs” directory stores all the files generated when the pipeline is used. 

It includes the final output files:
- **summary_output.csv**
- **detailed_output.csv** 
- **data_provenance.txt** 

and intermediate files, such as the ones with details on the knock-ins alternatives for growth on uncommon substrates (**consumption.csv**) and the .tsv files for MDF calculation. 

### summary_output.csv

**summary_output.csv** contains the score of each criterion per model, along with the final score, the ranking and the percentage of matching Biobrick part. If Optknock was performed, the information on the reaction knock-outs is also included in this file. 

[!summary ouptut](/outputs/summary.png)

The figure above shows an example output. In the large table, the first five columns identify the model of the knock-in variants and report their associated scores per criteria. The next column indicates the percentage of matching parts per model variant. The last three column report the result of Optknock. The second table below the one just mentioned shows how the model variants are ranked according to the evaluation of their properties with the introduced knock-ins.

### detailed_output.csv

**detailed_output.csv** specifies the general settings of the analysis (i.e. the constraints used in the sequential optimizations, the equation of the conversion and the expression host). It also includes all the information on the fluxes of uptake and production of each alternative knock-in strategy, the pathway’s MDF value, the identifiers of the heterologous reactions and the information on the associated Biobricks, including the harmonised sequences.

[!detailed output](/outputs/detailed.png)

The figure above shows an example of the detailed output file.The smaller box contains information on the general settings of the analysis. These includes the constraints used during the optimization of growth on the selected substrate (optimization of the uptake rate) and the ones used during the optimization of production rate. In addition the balanced equation and the species and model of the organism of choice are indicated as written in the input file. The larger table contains the detained information specific per engineering strategy: fluxes of consumption and production, values per every criteria, EC numbers of the heterologous reactions, their associated parts if found, the native and RSs-free harmonised sequences.