#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Snakemake template python snakemkake script for testing

"""
import logging
from utils import mock_snakemake
import pandas as pd


if __name__ == "__main__":

    if 'snakemake' not in globals():

        snakemake = mock_snakemake(
            rulename='test_rule',  # change if rule is named differently
            case="1",  # change to required wildcards
        )

    #helper.configure_logging(snakemake)

    logger = logging.getLogger(__name__)
    
    #%% load data
    #config = helper.load_config_yaml(snakemake.input.config)
    test_data = pd.read_csv(snakemake.input.input_file)
    
    print(f'Writing file {snakemake.output.output_file}')
    with open(snakemake.output.output_file,'w') as fid:
        fid.writelines(['test'])