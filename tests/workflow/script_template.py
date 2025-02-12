#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Snakemake template python snakemkake script for testing

"""
import logging
from utils import mock_snakemake
import pandas as pd
import os

if __name__ == "__main__":

    if 'snakemake' not in globals():

        snakemake = mock_snakemake(
            rulename='rule_auto_created',  #SNAKE_DEF
            case="1",  # change to required wildcards
        )

    #SNAKE_DEF: environment :  "envs/environment.yml"

    logger = logging.getLogger(__name__)
    
    #%% load data
   
    config = os.readlines(snakemake.input.config) #SNAKE_DEF: 'configfile'
    test_data = pd.read_csv(snakemake.input.input_file_1a) #SNAKE_DEF: 'input_{case}.csv'
    
    #%% calculation

    
    #%% writing output
    
    print(f'Writing file {snakemake.output.output_file}')  
    with open(snakemake.output.output_file,'w') as fid: #SNAKE_DEF: 'output_{cast}.csv'
        fid.writelines(['test'])
