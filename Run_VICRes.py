#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 11:13:54 2024

@author: dardiry
"""

import math
import pandas as pd
import subprocess
import os
import numpy as np
import calendar
import glob

# %% Define Working Directory
wd='/Users/dardiry/Academia/Cornell/Research/Integrated_Modeling/Model_Couplers/Analysis/VIC/Software_Package/VICRes_StepByStep_Hisham/Toy_Example/Routing/SourceCode/'
os.chdir(wd)

results_dir='/Users/dardiry/Academia/Cornell/Research/Integrated_Modeling/Model_Couplers/Analysis/VIC/Software_Package/VICRes_StepByStep_Hisham/Toy_Example/Results/'
# %% Read Configuration File and Update the Simulation Step

config_file ='../../RoutingSetup/configuration.txt'
fortran_executable='./rout' 


with open(config_file, 'r') as file:
        lines = file.readlines()
        
stepbystep_line = 33  # line number in the configuration file
lines[stepbystep_line-1] = '.true.' + '\n'

sim_day_line = 35  # line number in the configuration file

START_YEAR=2000
END_YEAR=2000
START_MONTH=1
END_MONTH=4

k=0

for year in range(START_YEAR,END_YEAR+1):
    for month in range(START_MONTH,END_MONTH+1):
        num_days=calendar.monthrange(year, month)
        for day in range(1,num_days[1]+1):
            
            new_step = '%d %d %d %d %d %d'%(START_YEAR,START_MONTH,1,year,month,day)
            lines[sim_day_line-1] = new_step
            # Write the modified lines back to the file
            with open(config_file, 'w') as file:
                file.writelines(lines)
            k=k+1
            print('Running VIC-Res for Step %d: %04d-%02d-%02d'%(k,year,month,day))    
            try:
                result = subprocess.run([fortran_executable,config_file], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                print(f"Program output:\n{result.stdout}")
            except subprocess.CalledProcessError as e:
                print(f"Error running Fortran program:\n{e.stderr}")    
                
# %% Remove Intermediate Files
directory=results_dir+'STEPBYSTEP/'
pattern='*STEP*'
files_to_remove = glob.glob(os.path.join(directory, pattern))
for file_path in files_to_remove:
    try:
        os.remove(file_path)
        print(f'{file_path} has been removed successfully.') 
    except FileNotFoundError:
        print(f'{file_path} does not exist.')
    except PermissionError:
        print(f'Permission denied: {file_path}')
    except Exception as e:
        print(f'Error removing {file_path}: {e}')



