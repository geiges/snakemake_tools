
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimal template of scriptfile

"""

#%% Imports
from snakemake_tools.helper import mock_snakemake, configure_logging, update_snakefile
import logging
#%% Function defintion
def do_something(input_data):
    pass
    


#%% Main execution

if __name__ == '__main__':
    #%% config
    
    

    if 'snakemake' not in globals():
        
        # optional -> comment in if snakefile should be updated each time this file is executed
        #update_snakefile('snakefile',  __file__)
    
        snakemake = mock_snakemake(
            rulename='rule_name',  #SNAKE_DEF
            wildcard_1="wc1",  # change to required wildcards
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
    
    
