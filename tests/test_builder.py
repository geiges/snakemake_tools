#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 09:01:04 2025

@author: andreasgeiges
"""

import os
def test_standart_snakemake():
    
    try:
        os.remove('result1.dummy')
    except OSError:
        pass


    os.system('snakemake -c1')    
    
    assert os.path.exists('result1.dummy')
    os.remove('result1.dummy')
    
def test_add_rule():
    

    