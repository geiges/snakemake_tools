# snakemake_tools
Tools to work more efficient with Snakemake

This package aims to simplify the use of snakemake in combination with python 
script files. 

## Snakefile builder
The builder tools allows to define all components of a snakemake rule within the
script file and provides a routine that automatically generates the rule blocks in
the snakefiles or allows to automatically update those.

The tag #SNAKE_DEF is used in all relevant lines that are necessary to create the rules
automatically.

The following minimal example containes all necessary information as an example to generate
the rule in the snakefile by 
snakemake_tools.builder.add_rule_to_snakefile(
    snakefile=
    script_file=
    )
    
Complex inputs using functions or similar things require manual edits still, while this
automatisation only covers the basic workfow and avoids tedious manual copy/pasting.




''' Minimal script template
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimal template of scriptfile

"""
from snakemake_tools.helper import mock_snakemake, configure_logging

# Function defintion
def do_something(input_data):
    pass
    
# Main execution

if __name__ == '__main__':
    #%% config
    
    

    if 'snakemake' not in globals():
        
        update_snakefile('snakefile',  __file__)
    
        snakemake = mock_snakemake(
            rulename='climate_assessment',  #SNAKE_DEF
            case="basecase",  # change to required wildcards#
            infilling='eqw',
        )
        
         #SNAKE_DEF logfile : '../output/{case}_{infilling}/log_climate_assessment.log'
         #SNAKE_DEF environment :  "../envs/environment.yml"
         
         logger = logging.getLogger(__name__)
         configure_logging(snakemake)
         
         input_data_path = snakemake.input.input_var_name #SNAKE_DEF 'cases/{case}/path/to/input.data'
        
        
    #%% Computation
    
    do_something(input_data_path)
    
    #%% Output saving
    output_data_path  = snakemake.output.output_var_name #SNAKE_DEF 'cases/{case}_path_to_output.data'
    
    
'''

## Snakemake helper
The helper tools provides some convenience tools to work with python. 


# Contribution
Some of the helper functions are adaped from the package "pypsa_eur"  (see https://github.com/PyPSA/pypsa-eur) and are based on their 
incredible work.