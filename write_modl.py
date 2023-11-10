"""" """


import numpy as np
import pandas as pd
import os
import shutil
# import sys
import re
import tempConstants as tc
# sys.path.append("C:\\Users\\jgazquez\\RealSONIC")
import tempFunctions as tf

""""--------------- INPUTS ---------------"""
cell_nr = 7
cell_folder = "L23_PC_cADpyr229_2"
sec_type = 'somatic'
pickle_folder = "c:\\users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\test_joa\\"
pickle_file = "realneuron_lookups_fs1.00.pkl"

""""--------------------------------------"""

mech_folder = "cells/"+cell_folder+"/mechanisms/" #"mechanisms\\"
mod_files, mod_names = tf.read_mod(mech_folder) #returns a list with mod file content and the names of the mod files
l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names) #one line formulas but these aren't used, only the 'hit' lines
states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless, gives a dictionary of all the state names, as can be seen in PySONIC\neurons\real_neuron.py

#pkl_txt = tf.read_pickle(pickle_folder,pickle_file)
pkl_txt = pd.read_pickle(pickle_folder+pickle_file) #read pickle file
func_tables = pkl_txt['tables'].keys() #this returns all the gating state kinetics, as can be seen in PySONIC\neurons\real_neuron.py
# for e in func_tables:
#     print(e)

mod_files = []
mod_names = []
for root, dirs, files in os.walk(mech_folder): #go through all files in the mechanics folders (all depths)
    for file in files:
        if file.endswith(".mod") and not 'eff' in root+file: #we only want to duplicate mechanism .modl files and not already created effective copies
            file_dupl = file.replace(".mod","_eff.mod") #we create the effective 'duplicate'which is currently empty
            # first we copy everything from .mod to _eff.mod without the PROCEDURE rates() block
            block = None #block keeps track in which BLOCK the writer is at the moment
            mod_eff = False #only remove and add blocks if the mechanic is voltage dependent, otherwise just copy
            for e in func_tables:
                if file.replace(".mod","").replace("_","") in e:
                    mod_eff = True #the mechanic is voltage dependent and needs an effective duplicate
            if not mod_eff: #just copy the file if the mechanism is not voltage dependent and go to the next mechanism (hence the continue)
                shutil.copy(root+file,root+"eff\\"+file) #remove the _eff when no effective variables are pretabulated to indicate that it is a pure duplicate without adaptations
                continue

            #start writing to the new file
            with open (root+file) as f, open(root+"eff\\"+file_dupl,'w') as dupl: #now the effective (still empty) duplicate will copy everything except the PROCEDURE block
                flist = list(f)
                for i, line in enumerate(flist):
                    if re.search('v.*\(m[Vv]\)',line): #line needs to be replaced
                        dupl.write('\tv (nC/cm2)\n\tVm (mV)\n') #add/replace specific lines -> recast
                        continue
                    elif re.search('i.*\=.*v',line): #line needs to be replaced
                        voltages_hits = re.findall('\Wv\W',line) #\W means: "no word" do v cannot be preceded or followed by a combination of letters
                        if len(voltages_hits) > 1:
                            print(f'{tc.bcolors.OKCYAN} replaced more than 1 time v by Vm in {file}\n -> {line} {tc.bcolors.ENDC}')
                        for e in voltages_hits:
                            line = line.replace(e,e.replace('v','Vm')) #this is the actual potential and not the v that is used in NEURON
                        dupl.write(line) #add specific lines -> recast
                        continue
                    elif block == "INITIAL" and '=' in line: #line/equations needs to be replaced
                        LHS,RHS = line.split('=') #split equation in LHS/RHS
                        var = re.findall('[a-zA-Z]',LHS)[0]
                        alph = f"alpha{var}_{file.replace('.mod','').replace('_','')}(A_t, y)"
                        bet = f"beta{var}_{file.replace('.mod','').replace('_','')}(A_t, y)"
                        dupl.write(f"{LHS}= {alph} / {alph} + {bet}\n") #all gating variables have the same type of formula # see PySONIC/neurons/real_neurons.py in steadyStates
                        continue
                    elif block == "DERIVATIVE" and '=' in line: #line/equations needs to be replaced
                        LHS,RHS = line.split('=') #split equation in LHS/RHS
                        var = re.findall('[a-zA-Z]\'',LHS)[0][:-1] #here we remove the differential ' from the actual variable
                        alph = f"alpha{var}_{file.replace('.mod','').replace('_','')}(A_t, y)"
                        bet = f"beta{var}_{file.replace('.mod','').replace('_','')}(A_t, y)"
                        dupl.write(f"{LHS}= {alph} * (1 - {var}) - {bet} * {var}\n") #all gating variables have the same type of formula # see PySONIC/neurons/real_neurons.py in derStates
                        continue
                    if re.search('^[A-Z][A-Z]*.*\{',line): #do determine if we are in a specific block or not
                        block = re.search('^[A-Z][A-Z]*',(re.search('^[A-Z][A-Z]*.*\{',line).group(0))).group(0) #first look if we are in a BLOCK initiation line and then extract the actual block
                        #print(block)
                    if block == "PROCEDURE" or re.search('rates()',line): #do not copy the PROCEDURE block or the rates() function caller -> is replaced with FUNCTION TABLES and called as can be seen above in the INITIAL and DERIVATIVE blocks
                        continue
                    dupl.write(line) #copy line if all cases above are not the case
                
                    if block == "NEURON" and flist[i+1].startswith('}'): #put extra lines at the end of the NEURON block
                        dupl.write("\tRANGE Adrive, Vm, y, Fdrive, A_t : section specific\n\tRANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)\n") #add specific lines
                    elif block == "PARAMETER" and flist[i].startswith("PARAMETER"): #put extra lines at the beginning of the PARAMETER block
                        dupl.write("\tstimon       : Stimulation state\n\tFdrive (kHz) : Stimulation frequency\n\tAdrive (kPa) : Stimulation amplitude\n\tdetailed     : Simulation type\n") #add specific lines
                    elif block == "ASSIGNED" and flist[i+1].startswith('}'): #add extra lines at the end of the ASSIGNED block
                        dupl.write("\tA_t  (kPa)\n\ty\n") #add specific lines
                    elif (block == "BREAKPOINT" and flist[i].startswith("BREAKPOINT")) or (block == "INITIAL" and flist[i].startswith("INITIAL")): #add extra line at the beginning of these 2 blocks
                        dupl.write("\tupdate()\n") #add specific line
                        
                    if line.startswith('}'):
                        if block == "ASSIGNED": #put the include and the FUNCTION TABLES after the assigned
                            dupl.write("\nINCLUDE \"update.inc\"\n\n")  #include this file
                            for e in func_tables: #check if it is an effective 'duplicate'
                                if len(e) <= 2 or file.replace('.mod','').replace('_','') in e: # 'if len(e) <= 2' is used to always include V (which is only 1 char)
                                                                                                # if the mechanism is in the gating state kinetic, include it also as a FUNCTION TABLE
                                    # then we append the FUNCTION TABLE lines/blocks
                                    dupl.write("FUNCTION_TABLE ")
                                    dupl.write(e)
                                    dupl.write("(A(kPa), Q(nC/cm2)) (mV)\n")    
                        block = None


# from neuron import h
# h.load_file('init.hoc')
# h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
# h.cell_chooser(cell_nr)
