#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 09:01:04 2025

@author: andreasgeiges
"""
import pathlib
import shutil 
import os
from snakemake_tools import builder
#%%
TEST_PATH = pathlib.Path(__file__).parent 

SNAKEFILE_TEMPLATE = TEST_PATH  /  'workflow/snakefile'
SNAKEFILE_TEMP = '/tmp/snakefile'
def test_standart_snakemake():
    
    try:
        os.remove('result1.dummy')
    except OSError:
        pass

    os.chdir(TEST_PATH)
    os.system('snakemake -c1')    
    
    assert os.path.exists('result1.dummy')
    os.remove('result1.dummy')
    
def test_add_rule():
    shutil.copy(SNAKEFILE_TEMPLATE, SNAKEFILE_TEMP)
    builder.add_snakemake_rule(snakefile=SNAKEFILE_TEMP, 
                              scriptfile= TEST_PATH  / 'workflow/script_template.py')
    
    with open(SNAKEFILE_TEMP,'r') as fid:
        content_obs = fid.read()
        
    with open( TEST_PATH  / 'workflow/snakefile_exp_test_add','r') as fid:
        content_exp = fid.read()
        
    assert content_exp == content_obs

def test_update_rule():
    
    shutil.copy(SNAKEFILE_TEMPLATE, SNAKEFILE_TEMP)
    builder.add_snakemake_rule(snakefile=SNAKEFILE_TEMP, 
                              scriptfile= TEST_PATH  / 'workflow/script_template.py')
    
    _ = builder.update_snakefile(snakefile=SNAKEFILE_TEMP, 
                              scriptfile= TEST_PATH  / 'workflow/script_template_updated.py')

    with open(SNAKEFILE_TEMP,'r') as fid:
        content_obs = fid.read()
        
    with open( TEST_PATH  /'workflow/snakefile_exp_test_update','r') as fid:
        content_exp = fid.read()
        
    assert content_exp == content_obs

#test_update_rule()