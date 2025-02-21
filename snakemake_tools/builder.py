#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 08:58:00 2025

@author: andreasgeiges
"""
import re
import pathlib
import datetime
import shutil


def _find_rulename(string):
    
    if '\'' in string:
        rulename = re.findall("'(.*?)'", string)[0]
    elif '\"' in string:
        rulename = re.search('"(.*?)"', string)[0]
    else:
        rulename = None
    return rulename

def _find_input_def(string):
    
    #if string.endswith("'\n")
    if string.endswith("'\n") and ('#SNAKE_DEF' in string):
        input_file = re.findall("#SNAKE_DEF.*'(.*?)'", string)[0]
        input_name = re.findall("snakemake.input.([a-zA-Z_$0-9]*)", string)[0]
        input_file = f"'{input_file}'"
    elif string.endswith('"\n'):
        input_file = re.findall('#SNAKE_DEF.*"(.*?)"', string)[0]
        input_name = re.findall("snakemake.input.([a-zA-Z_$0-9]*)", string)[0]
        input_file = f"'{input_file}'"
    elif string.endswith(")\n") or string.endswith("]\n"):
        #assuming functional input
        
        input_file = re.findall("#SNAKE_DEF (.*?)\n", string)[0]
        input_name = re.findall("snakemake.input.([a-zA-Z_$0-9]*)", string)[0]

    else:
        pass
        #input_file = None
    return input_file, input_name

def _find_output_def(string):
    
    if string.endswith("'\n") and ('#SNAKE_DEF' in string):
        output_file = re.findall("#SNAKE_DEF.*'(.*?)'", string)[0]
        output_file = f"'{output_file}'"
        
    elif string.endswith('"\n'):
        output_file = re.findall('#SNAKE_DEF.*"(.*?)"', string)[0]
        output_file = f"'{output_file}'"
    
    elif string.endswith(")\n") or string.endswith("]\n"):
        #assuming functional input
        output_file = re.findall("#SNAKE_DEF (.*?)\n", string)[0]
        output_file = f"'{output_file}'"
    else:
        pass
        #output_file = None
    output_name = re.findall("snakemake.output.([a-zA-Z_$0-9]*)", string)[0]
    return output_file, output_name

def find_environment(string):
    if ('\'' in string):
        match_str = re.findall("#SNAKE_DEF.*'(.*?)'", string)[0]
        # output_name = re.findall("snakemake.output.(.*?)\)", string)[0]
    elif '\"' in string:
        match_str = re.findall('#SNAKE_DEF.*"(.*?)"', string)[0]
        # output_name = re.findall("snakemake.output.(.*?)\)", string)[0]
    else:
        match_str =  None
    return match_str

def find_log(string):
    if ('\'' in string):
        match_str = re.findall("#SNAKE_DEF.*'(.*?)'", string)[0]
        # output_name = re.findall("snakemake.output.(.*?)\)", string)[0]
    elif '\"' in string:
        match_str = re.findall('#SNAKE_DEF.*"(.*?)"', string)[0]
        # output_name = re.findall("snakemake.output.(.*?)\)", string)[0]
    else:
        match_str =  None
    return match_str
    
def  get_rule_string_block(scriptfile):
     input_files = list()
     output_files = list()
     environment = None
     log_file = None
     with open(scriptfile, 'r') as fid:
         
         for line in fid.readlines():
             #print(line)
             
             if '#SNAKE_DEF' in line:
             
                 if 'rulename' in line:
                     
                     rulename = _find_rulename(line )
                     
                     
                 if  'snakemake.input.' in line:
                     input_file, input_name = _find_input_def(line )
                     if input_file is not None:
                         input_files.append(f'{input_name} = {input_file},')
                     
                 if  'snakemake.output.' in line:
                     output_file, output_name = _find_output_def(line )
                     
                     if output_file is not None: 
                         output_files.append(f'{output_name} = {output_file},') 
                 if "environment" in line:
                     environment = find_environment(line)
                     
                 if "logfile" in line:
                     log_file = find_log(line)
                         
                 
         #print(input_files)
         #print(output_files)
         
         # Concatenate rule string
         rule_strings = list()
         rule_strings.append(f'rule {rulename}:')
        
         rule_strings.append('\tinput:')
         for inputstr in input_files:
             rule_strings.append(f'\t\t{inputstr}')
         
         if environment is not None:
             rule_strings.append('\tconda:')
             rule_strings.append(f'\t\t"{environment}"')
         
         if log_file is not None:
             rule_strings.append('\tlog:')
             rule_strings.append(f'\t\tlog = "{log_file}"')
         
         rule_strings.append('\toutput:')
         for outputstr in output_files:
             rule_strings.append(f'\t\t{outputstr}')
         rule_strings.append('\tscript:')
         rule_strings.append(f'\t\t"{pathlib.Path(scriptfile).name}"')
        
         rule_strings.append('\n')
         
     return rule_strings
     
def add_snakemake_rule(snakefile, scriptfile):
    
    
    rule_strings = get_rule_string_block(scriptfile)
    
    with open(snakefile, "a") as myfile:
        myfile.write('\n'.join(rule_strings))

         
def update_snakefile(snakefile,
                     scriptfile):
    
    filename = pathlib.Path(scriptfile).name
    print('Making backup of snakefile to /tmp')
    shutil.copyfile(snakefile, f'/tmp/{filename}_{datetime.datetime.now().strftime("%Y%m%d_%H%m%S")}.bkp')
    with open(snakefile, 'r') as fid:
        org_lines = fid.readlines()
    
    rule_strings = get_rule_string_block(scriptfile)
    print(rule_strings)
    with open(snakefile,'w') as fid:
        
        skip_line=False
        found = False
        for line in org_lines:
            
            
            if  line.startswith('rule') and skip_line:
                #stop  skipping org input
                skip_line = False
                
                
                
                
            if rule_strings[0] in line:
                
                fid.write('\n'.join(rule_strings))
                #start skipping org input
                skip_line =True
                found = True
                
            if not skip_line:
                fid.write(line)
                
    if not found:
        print('rule not found, adding new rule')
        add_snakemake_rule(snakefile, scriptfile)
            
           
            
            
        
    

# string = update_snakefile(snakefile='../tests/workflow/snakefile', 
#                           scriptfile='../tests/workflow/script_template.py')


