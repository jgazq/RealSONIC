"""" """


import numpy as np
import pandas as pd
import os
import shutil
import sys
sys.path.append("C:\\Users\\jgazquez\\RealSONIC")
import tempFunctions as tf

""""--------------- INPUTS ---------------"""
cell_nr = 7
cell_folder = "L23_PC_cADpyr229_2"
sec_type = 'somatic'
pickle_folder = "c:\\users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\test_joa\\"
pickle_file = "realneuron_lookups_fs1.00.pkl"

""""--------------------------------------"""

mech_folder = "cells/"+cell_folder+"/mechanisms/" #"mechanisms\\"
mod_files, mod_names = tf.read_mod(mech_folder)
l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names)
states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless

#pkl_txt = tf.read_pickle(pickle_folder,pickle_file)
pkl_txt = pd.read_pickle(pickle_folder+pickle_file)
func_tables = pkl_txt['tables'].keys()
for e in func_tables:
    print(e)

mod_files = []
mod_names = []
for root, dirs, files in os.walk(mech_folder):
    for file in files:
        if file.endswith(".mod") and not 'eff' in root+file:
            file_dupl = file.replace(".mod","_eff.mod")
            # first we copy everything from .mod to _eff.mod without the PROCEDURE rates() block
            copy = True
            mod_eff = False #only remove and add blocks if the mechanic is voltage dependent, otherwise just copy
            for e in func_tables:
                if file.replace(".mod","") in e:
                    mod_eff = True
            if not mod_eff: #just copy the file
                print(file)
                shutil.copy(root+file,root+"eff\\"+file) #remove the _eff when no effective variables are pretabulated
                continue
            with open (root+file) as f, open(root+"eff\\"+file_dupl,'w') as dupl:
                for line in f:
                    if line.startswith('PROCEDURE'):
                        copy = False
                    if copy:
                        dupl.write(line)
                    if line.startswith('}'):
                        copy = True
            # then we append the FUNCTION TABLE lines/blocks
            with open (root+file) as f, open(root+"eff\\"+file_dupl,'a') as dupl:
                for e in func_tables:
                    if len(e) <= 2 or file.replace('.mod','') in e:
                        dupl.write("FUNCTION_TABLE ")
                        dupl.write(e)
                        dupl.write("(A(kPa), Q(nC/cm2)) (mV)\n")           
                


# from neuron import h
# h.load_file('init.hoc')
# h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
# h.cell_chooser(cell_nr)
