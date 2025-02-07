"""" this files rewrites all the mechanism .mod files in a mechanisms folder to replace the gating kinetics functions with the lookup tables """


import numpy as np
import pandas as pd
import os
import shutil
# import sys
import re
import tempConstants as tc
# sys.path.append("C:\\Users\\jgazquez\\RealSONIC")
import tempFunctions as tf
from neuron import h 
h.load_file("init.hoc")

""""--------------- INPUTS ---------------"""
cell_nr = 7
sec_type = 'somatic'
pickle_folder = "/Users/jgazquez/PySONIC/PySONIC/lookups/test_joa/1overtone/"#"c:\\users\\jgazquez\\PySONIC\\PySONIC\\lookups\\"
pickle_file = "realneuron_lookups_1overtone_Jac.pkl" #overtones = 2 will only work if this is also lookups_2overtone.pkl
DEBUG = 0
Cm0_var = 1
overtones = 1
VERBATIM = 1

""""--------------------------------------"""

"""load in state names from the mechanisms folder and their gating state kinetics from the .pkl LUT"""
cell_folder = h.cell_names[cell_nr-1].s
mech_folder = "mechanisms/"#"cells/"+cell_folder+"/mechanisms/" #"mechanisms\\"
mod_files, mod_names = tf.read_mod(mech_folder) #returns a list with mod file content and the names of the mod files
l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names) #one line formulas but these aren't used, only the 'hit' lines
states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless, gives a dictionary of all the state names, as can be seen in PySONIC\neurons\real_neuron.py

#pkl_txt = tf.read_pickle(pickle_folder,pickle_file)
pkl_txt = pd.read_pickle(pickle_folder+pickle_file) #read pickle file
func_tables = pkl_txt['tables'].keys() #this returns all the gating state kinetics, as can be seen in PySONIC\neurons\real_neuron.py
# print(func_tables)
# [print(e) for e in func_tables]
v_str = 'vbt' if VERBATIM else ''

# tf.mod_duplicate(mech_folder,[e.split("_")[-1] for e in func_tables]);quit() # _eff.mod are copies from the original ones -> this is only used to see git changes

"""iterate over all mechanisms files and create duplicates"""
mod_files = []
mod_names = []
for root, dirs, files in os.walk(mech_folder): #go through all files in the mechanics folders (all depths)
    for file in files:
        if file.endswith(".mod") and not 'eff' in root+file: #we only want to duplicate mechanism .modl files and not already created effective copies
            file = tf.one_to_multiline(root,file) #replace all the BLOCKS that are defined in one line to a block that is defined on multiple lines in order to comply with future adaptations
            file_repl = file.replace(".mod","").replace("_","") #a pruned version of the mod name used for recognition later

            #overtone variables
            overtone_NEURON = ''
            overtone_ASSIGNED = ''
            overtone_FUNCTION_TABLE = ''
            A_LUT = ['V']
            B_LUT = []
            overtone_ARGUMENTS = ''

            for overtone in range(overtones):
                overtone_NEURON += f', a{overtone+1}, b{overtone+1}' #add 'qi, fi' for every overtone
                overtone_ASSIGNED += f'\ta{overtone+1}  (nC/cm2)\n'
                overtone_ASSIGNED += f'\tb{overtone+1}  (rad)\n'

            # first we copy everything from .mod to _eff.mod without the PROCEDURE rates() block
            block = None #block keeps track in which BLOCK the writer is at the moment
            mod_eff = False #only remove and add blocks if the mechanic is voltage dependent, otherwise just copy
            for e in func_tables:
                if e.endswith(file_repl):
                    mod_eff = True #the mechanic is voltage dependent and needs an effective duplicate
            file_dupl = file.replace(".mod",f"_eff_{overtones}ov.mod") if mod_eff else file #we create the effective 'duplicate'which is currently empty

            """"the mechanisms, where there are no computed effective variables and function tables, also need some modifications for the recasting"""
            """except xtra.mod?"""
            if 'xtra' in file_repl or 'CaDynamics' in file_repl:
                shutil.copy(root+file,root+f"eff_{overtones}ov{v_str}\\"+file) #remove the _eff when no effective variables are pretabulated to indicate that it is a pure duplicate without adaptations #DONT DO THIS IN MAC
                continue    
            if Cm0_var:
                if 'pas' in file_repl:
                    #here we don't iterate over the Cm0_map as it doesn't include the value for 0.02 (which is not in LUT)
                    "Cm0 = 1"
                    #shutil.copy(root+file,root+f"eff_{overtones}ov\\"+file.replace(".mod",f"_eff_{overtones}ov.mod")) #DONT DO THIS IN MAC
                    with open(os.path.join(root,file)) as f, open(os.path.join(root,f"eff_{overtones}ov{v_str}",file.replace(".mod",f"_eff_{overtones}ov.mod")),'w') as dupl: 
                        flist = list(f)
                        flist_ov = tf.add_custom_pas(flist,overtones,1,vbt=VERBATIM)
                        dupl.writelines(flist_ov)
                    "Cm0 = 2"
                    with open(os.path.join(root,file)) as f, open(os.path.join(root,f"eff_{overtones}ov{v_str}",file.replace(".mod",f"_eff_{overtones}ov.mod").replace('.mod','_2.mod')),'w') as dupl: 
                        flist = list(f)
                        flist2 = tf.SUFFIX_Cm0(flist,"2")
                        flist2_ov = tf.add_custom_pas(flist2,overtones,2,vbt=VERBATIM)
                        dupl.writelines(flist2_ov)
                    "Cm0 = 0.02"
                    with open(os.path.join(root,file)) as f, open(os.path.join(root,f"eff_{overtones}ov{v_str}",file.replace('.mod','_0_02.mod')),'w') as dupl: 
                        flist = list(f)
                        flist = tf.eff_to_noteff(flist,0.02)
                        flist02 = tf.SUFFIX_Cm0(flist,"0_02")
                        flist02_ov = tf.add_custom_pas(flist02,overtones,0.02,vbt=VERBATIM)
                        dupl.writelines(flist02)
                    continue   


            # if not mod_eff: #just copy the file if the mechanism is not voltage dependent and go to the next mechanism (hence the continue)
            #     shutil.copy(root+file,root+"eff\\"+file) #remove the _eff when no effective variables are pretabulated to indicate that it is a pure duplicate without adaptations
            #     continue      
  
            """start writing to the new file"""
            with open(os.path.join(root,file)) as f, open(os.path.join(root,f"eff_{overtones}ov{v_str}",file_dupl),'w') as dupl: #now the effective (still empty) duplicate will copy everything except the PROCEDURE block
                flist = list(f)
                voltage_gated = False #we wait for proof to assume the mechanism is voltage gated
                no_indep = True #until we encounter an INDEPENDENT block, there is none
                
                for i, line in enumerate(flist):
                    if re.search(tc.onelineBLOCK_pattern,line):
                        print(f"{line} in {file}")
                    if block == "PROCEDURE": #do not copy the PROCEDURE block or the rates() function caller -> is replaced with FUNCTION TABLES and called as can be seen above in the INITIAL and DERIVATIVE blocks
                        if mod_eff:
                            continue
                        else:
                            dupl.write(line)
                            continue
                    elif re.search(tc.neuronvoltage_pattern,line) and re.findall(tc.vsecluded_pattern,line): #line needs to be replaced #for both noneff and eff
                        dupl.write('\tv (nC/cm2)\n\tVm (mV)\n') #add/replace specific lines -> recast
                        voltage_gated = True #if mechanism contains v, it is assumed to be voltage gated and needs to be recast
                        continue
                    elif re.search(tc.vRHS_pattern,line) and re.findall(tc.vsecluded_pattern,line):#elif re.search('i.*\=.*v',line): #line needs to be replaced #for both
                        voltages_hits = re.findall(tc.vsecluded_pattern,line) #\W means: "no word" do v cannot be preceded or followed by a combination of letters
                        if len(voltages_hits) > 1:
                            print(f'{tc.bcolors.OKCYAN} replaced more than 1 time v by Vm in {file}\n -> {line} {tc.bcolors.ENDC}')
                        for e in voltages_hits:
                            line = line.replace(e,e.replace('v','Vm')) #this is the actual potential and not the v that is used in NEURON
                        dupl.write(line) #add specific lines -> recast
                        continue
                    elif (block == "INITIAL" and '=' in line) and mod_eff: #line/equations needs to be replaced
                        LHS,RHS = line.split('=') #split equation in LHS/RHS
                        var = re.findall(tc.state_pattern,LHS)[0]
                        alph = f"alpha{var}_{file_repl}(A_t, y{overtone_ARGUMENTS})"
                        bet = f"beta{var}_{file_repl}(A_t, y{overtone_ARGUMENTS})"
                        if DEBUG:
                            #dupl.write(f'printf("{file}: V = %g, alpha = %g, beta = %g\\n",V(A_t,y), {alph}, {bet})\n') #add line for debugging
                            dupl.write(f'printf("{file}: \\n")\n') #add line for debugging
                            dupl.write(f'printf("V = %g\\t",V(A_t,y))\n')
                            dupl.write(f'printf("alpha = %g\\t" ,{alph})\n')
                            dupl.write(f'printf("beta = %g\\n" ,{bet})\n')
                            #dupl.write(f'printf("V = %g, alpha = %g, beta = %g\\n",V(A_t,y), {alph}, {bet})\n') #add line for debugging
                        dupl.write(f"{LHS}= {alph} / ({alph} + {bet})\n") #all gating variables have the same type of formula # see PySONIC/neurons/real_neurons.py in steadyStates
                        continue
                    elif (block == "DERIVATIVE" and '=' in line) and mod_eff: #line/equations needs to be replaced
                        LHS,RHS = line.split('=') #split equation in LHS/RHS
                        var = re.findall(tc.stateder_pattern,LHS)[0][:-1] #here we remove the differential ' from the actual variable
                        alph = f"alpha{var}_{file_repl}(A_t, y{overtone_ARGUMENTS})"
                        bet = f"beta{var}_{file_repl}(A_t, y{overtone_ARGUMENTS})"
                        dupl.write(f"{LHS}= {alph} * (1 - {var}) - {bet} * {var}\n") #all gating variables have the same type of formula # see PySONIC/neurons/real_neurons.py in derStates
                        continue
                    if re.search(tc.block_init_pattern,line): #do determine if we are in a specific block or not
                        block = re.search(tc.block_pattern,(re.search(tc.block_init_pattern,line).group(0))).group(0) #first look if we are in a BLOCK initiation line and then extract the actual block
                        #print(block)
                    if re.search('rates()',line) and mod_eff: #do not copy the PROCEDURE block or the rates() function caller -> is replaced with FUNCTION TABLES and called as can be seen above in the INITIAL and DERIVATIVE blocks
                        continue
                    dupl.write(line) #copy line if all cases above are not the case
                    if block == "NEURON" and flist[i+1].startswith('}'): #put extra lines at the end of the NEURON block
                        dupl.write(f"\tRANGE Adrive, Vm, y, Fdrive, A_t{overtone_NEURON} : section (even segment) specific\n\tRANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)\n") #add specific lines
                    if block == "UNITS" and flist[i+1].startswith('}'): #put extra lines at the end of the UNITS block
                        dupl.write("\tPI = (pi) (1) :in order to use the constant pi\n") #add specific lines
                    elif block == "PARAMETER" and flist[i].startswith("PARAMETER"): #put extra lines at the beginning of the PARAMETER block
                        dupl.write("\tstimon       : Stimulation state\n\tFdrive (kHz) : Stimulation frequency\n\tAdrive (kPa) : Stimulation amplitude\n\tdetailed     : Simulation type\n") #add specific lines
                    elif block == "ASSIGNED" and flist[i+1].startswith('}'): #add extra lines at the end of the ASSIGNED block
                        dupl.write(f"\tA_t  (kPa)\n\ty\n{overtone_ASSIGNED}") #add specific lines
                    elif ((block == "BREAKPOINT" and flist[i].startswith("BREAKPOINT")) or (block == "INITIAL" and flist[i].startswith("INITIAL"))) and voltage_gated: #and mod_eff: #add extra line at the beginning of these 2 blocks
                        if (block == "INITIAL" and flist[i].startswith("INITIAL")) and voltage_gated:
                            pass #in order to write v = init-value
                        if DEBUG:
                            if mod_eff:
                                #dupl.write(f'printf("{file}: V = %g\\n",V(A_t,y))\n') #add line for debugging
                                dupl.write(f'printf("{file}: \\n")\n') #add line for debugging
                                dupl.write(f'printf("V = %g\\n",V(A_t,y))\n') #add line for debugging
                        dupl.write("\tupdate()\n") #add specific line
                    elif (block == "INDEPENDENT"):
                        no_indep = False
                    
                    if line.startswith('}'):
                        if block == "ASSIGNED" and mod_eff and voltage_gated: #put the include and the FUNCTION TABLES after the assigned
                            dupl.write("\nINCLUDE \"update.inc\"\n")  #include this file
                            dupl.write("\n")
                            for e in func_tables: #check if it is an effective 'duplicate'
                                A = (e in A_LUT) #len(e) <= 2 #is it V or a LUT related to a variable of overtones
                                alphbet = e.endswith(file_repl)
                                B = (e in B_LUT)
                                if A or alphbet or B: # 'if len(e) <= 2' is used to always include V (which is only 1 char)
                                                                                                # if the mechanism is in the gating state kinetic, include it also as a FUNCTION TABLE
                                    # then we append the FUNCTION TABLE lines/blocks
                                    dupl.write("FUNCTION_TABLE ")
                                    dupl.write(e)
                                    dupl.write(f"(A(kPa), Q(nC/cm2){overtone_FUNCTION_TABLE}) (mV)\n") if (A or B) else dupl.write(f"(A(kPa), Q(nC/cm2){overtone_FUNCTION_TABLE}) (/ms)\n")
                        elif block == "ASSIGNED" and 'xtra' not in file_repl and voltage_gated:
                            dupl.write("\nINCLUDE \"update.inc\"") #_bis.inc\"\n")  #include this file
                            dupl.write("\n")
                            for e in A_LUT:
                                dupl.write("FUNCTION_TABLE ")
                                dupl.write(f"{e}(A(kPa), Q(nC/cm2){overtone_FUNCTION_TABLE}) (mV)\n") 
                            for e in B_LUT:
                                dupl.write("FUNCTION_TABLE ")
                                dupl.write(f"{e}(A(kPa), Q(nC/cm2){overtone_FUNCTION_TABLE}) (rad)\n") 
                        block = None
                if no_indep:
                    dupl.write('INDEPENDENT {\n\tt FROM 0 TO 1 WITH 1 (ms)\n}')
            """some lines need to be interchanged in the effective duplicate"""        
            with open(os.path.join(root,f"eff_{overtones}ov{v_str}",file_dupl),'r') as dupl:
                flist = list(dupl)
                new_list,worked = tf.str1_before_str2("LOCAL","update()",flist) #update() may not appear before LOCAL
                if worked:
                    flist = new_list
            with open(os.path.join(root,f"eff_{overtones}ov{v_str}",file_dupl),'w') as dupl:
                dupl.writelines(flist)
                    
                #print(any(['LOCAL' in e for e in flist]),file)
            if Cm0_var and voltage_gated:
                if "Prob" in file:
                    continue
                for Cm0fl, Cm0str in tc.Cm0_map.items():
                    if Cm0fl == 1: #do not touch file if Cm0 = 1
                        continue
                    with open(os.path.join(root,f"eff_{overtones}ov{v_str}",file_dupl),'r') as dupl:
                        flist = list(dupl)
                        "following lines: only to add a Cm0-suffix to the SUFFIX and the .mod debug prints"
                        #flistCm = tf.SUFFIX_Cm0(flist,Cm0str)
                        #flistCm = tf.replace_str(flist,'.mod',Cm0str+'.mod') #DEBUG line has to be specific for the Cm0 variant mechanism
                        "following lines: add a Cm0-suffix to all variables that contain the mechanism name"
                        torepl1 = file.replace('.mod','')
                        replwith1 = torepl1+Cm0str
                        torepl2, replwith2 = torepl1.replace('_',''), replwith1.replace('_','')
                        flistCm = tf.replace_str(flist,[torepl1, torepl2],[replwith1, replwith2]) 
                    with open(os.path.join(root,f"eff_{overtones}ov{v_str}",file_dupl.replace('.mod','_'+Cm0str+'.mod')),'w') as dupl_Cm0:
                        dupl_Cm0.writelines(flistCm)
print('Also update "update.inc" please!')
# from neuron import h
# h.load_file('init.hoc')
# h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
# h.cell_chooser(cell_nr)
