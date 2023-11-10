""""contains all the utility/auxiliary functions for both NEURONic as Pythonic"""

import matplotlib.pyplot as plt
# import PySONIC as ps
# import MorphoSONIC as ms
import numpy as np
# import pickle
import pandas as pd
from neuron import h
import re
import os
import copy
import shutil

import tempConstants as tc
import prev.Interp3Dfield as tt

# plot the distribution points of a transducer source
def plt_transdistr(psource,grid_type):
    nsources = psource.get_default_nsources()
    x, y = psource.getXYSources(m=nsources, d=grid_type)
    plt.title(f'{grid_type} distribution')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.plot(x * ms.M_TO_MM, y * ms.M_TO_MM, 'o', markersize = 4)
    plt.tight_layout()
    plt.show()

# read the LUT pickle files metadata
def read_pickle(path, filename):
    pkl = pd.read_pickle(path+filename)
    print('shape of V in tables',pkl['tables']['V'].shape)
    print('refs: ',pkl['refs'].keys())
    print('tables: ',pkl['tables'].keys())
    if 'tcomp' in pkl['tables'].keys(): #test files don't have a computation time value
        print('total simulation time: ',np.sum(pkl['tables']['tcomp'])/3600/24,'days')
    return pkl

# get the conductance bar values for every mechanism, the value is taken from the first encounter
def gbar_mechs():
    counter = [0,0,0,0,0,0,0]
    for sec in h.allsec():
        for seg in sec:
            for mech in seg:
                if str(mech) == "NaTa_t":
                    if counter[0] ==0:
                        print("NaTa_t")
                        print(mech.gNaTa_tbar)
                        counter[0] = 1
                if str(mech) == "Nap_Et2":
                    if counter[1] ==0:
                        print("Nap_Et2")
                        print(mech.gNap_Et2bar)
                        counter[1] =1
                if str(mech) == "NaTs2_t":
                    if counter[2] ==0:
                        print("NaTs2_t")
                        print(mech.gNaTs2_tbar)
                        counter[2] = 1
                if str(mech) == "K_Tst":
                    if counter[3] ==0:                
                        print("K_Tst")
                        print(mech.gK_Tstbar)
                        counter[3] = 1
                if str(mech) == "K_Pst":
                    if counter[4] ==0:   
                        print("K_Pst")
                        print(mech.gK_Pstbar)
                        counter[4] = 1
                if str(mech) == "Ih":
                    if counter[5] ==0:   
                        print("Ih")
                        print(mech.gIhbar)
                        counter[5] = 1
                if str(mech) == "Im":
                    if counter[6] ==0:
                        print("Im")
                        print(mech.gImbar)
                        counter[6] = 1

# print all mechanisms of a cell                        
def cell_mechs():
    mech_poss = []
    mech_actual = []

    mt = h.MechanismType(0)
    mname  = h.ref('')
    for i in range(mt.count()):
        mt.select(i)
        mt.selected(mname)
        mech_poss.append(mname[0])
    #mech_poss contains all possible mechanisms that this cell can have
    print(mech_poss)
    for e in mech_poss:
        for sec in h.allsec():
            if (h.ismembrane(e,sec=sec) and e not in mech_actual):
                print(e)
                mech_actual.append(e)
            for seg in sec:
                if (e not in mech_actual):
                    print(e)
                    mech_actual.append(e)   
                for mech in seg:
                    if (str(mech)==e and e not in mech_actual):
                        print("mech:  ",e)
                        mech_actual.append(e)                       
    print(mech_actual) 

# print every section with its corresponding mechanisms (every section contains segments
# that have the same mechanisms but are just spatially extended)                       
def sec_mechs():
    d_sec = {}
    prev = "9000"
    for sec in h.allsec():
        replace_sec = re.findall('\[[0-9]*\]',str(sec))
        sec_repl = str(sec)
        replace_prev = re.findall('\[[0-9]*\]',str(prev))
        prev_repl = str(prev)
        for e in replace_sec:
            sec_repl = sec_repl.replace(e,'')
        for e in replace_prev:
            prev_repl = prev_repl.replace(e,'')
        if sec_repl == prev_repl:
            prev = sec
            continue
        for seg in sec: 
            mech_list = []
            for mech in seg:
                mech_list.append(str(mech))
            sec_repl = 'somatic' if 'soma' in sec_repl else 'apical' if 'apic' in sec_repl\
                        else 'basal' if 'dend' in sec_repl else 'axonal' if 'axon' in sec_repl else sec_repl
            d_sec[sec_repl] = mech_list
            break
        prev = sec
    return d_sec

#prints the distance to the soma for all segments
def dist_sec_soma():
    i = 0
    for sec in h.allsec():
        for seg in sec:
            if "Scale" not in str(sec) and "Elec" not in str(sec):
                x_SOM,y_SOM,z_SOM = tt.findSOMA()
                dist_2_soma = np.sqrt((seg.x_xtra-x_SOM)**2 + (seg.y_xtra-y_SOM)**2 + (seg.z_xtra-z_SOM)**2)
                print(seg,dist_2_soma)

#returns a dict with all segments containing their 3D position
def coord_dict():
    i = 0
    dict = {}
    for sec in h.allsec():
        for seg in sec:
            if "Scale" not in str(sec) and "Elec" not in str(sec):
                dict[str(seg)] = (seg.x_xtra, seg.y_xtra, seg.z_xtra)
    return dict


#returns a dict with each type of segment containing all the segments and their locations
def coord_dict_type():
    i = 0
    dict = {}
    for sec in h.allsec():
        for seg in sec:
            replace_sec = re.findall('\[[0-9]*\]',str(sec))
            sec_repl = str(sec)
            for e in replace_sec:
                sec_repl = sec_repl.replace(e,'')
            sec_type = 'soma' if 'soma' in sec_repl else 'apical' if 'apic' in sec_repl\
            else 'basal' if 'dend' in sec_repl else 'axon' if 'axon' in sec_repl else sec_repl.lower()
            if "Scale" not in str(sec) and "Elec" not in str(sec):
                if sec_type not in dict.keys():
                    dict[str(sec_type)] = {str(seg) : (seg.x_xtra*1e-6, seg.y_xtra*1e-6, seg.z_xtra*1e-6)} #um to m
                else:
                    dict[str(sec_type)][str(seg)] = (seg.x_xtra*1e-6, seg.y_xtra*1e-6, seg.z_xtra*1e-6) #um to m           
    return dict

#read all .mod files in a directory and put them into a 2D-list
def read_mod(mech_folder, restrictions = None):
    mod_files = []
    mod_names = []
    for root, dirs, files in os.walk(mech_folder):
        for file in files:
            if restrictions:
                if file.endswith(".mod") and not 'eff' in root+file and file.replace(".mod","") in restrictions:
                    with open(root+file) as f:
                        lines_list = f.readlines()
                        mod_files.append(lines_list)
                        mod_names.append(file.replace('.mod','')) #root+file   
            elif not 'eff' in root+file:
                if file.endswith(".mod"):
                    with open(root+file) as f:
                        lines_list = f.readlines()
                        mod_files.append(lines_list)
                        mod_names.append(file.replace('.mod','')) #root+file   
    return mod_files, mod_names 

#read all .mod files in a directory and duplicate them in a seperate folder "eff"
def mod_duplicate(mech_folder, restrictions = None):
    mod_files = []
    mod_names = []
    for root, dirs, files in os.walk(mech_folder):
        for file in files:
            if restrictions:
                if file.endswith(".mod") and not file.endswith("_eff.mod") and file.replace(".mod","") in restrictions:
                    file_dupl = file.replace(".mod","_eff.mod")
                    shutil.copyfile(root+file,root+"eff\\"+file_dupl)
                    # with open(root+file_dupl,'w') as f:
                    #     lines_list = f.readlines()
                    #     mod_files.append(lines_list)
                    #     mod_names.append(file) #root+file   
            elif not file.endswith("_eff.mod"):
                if file.endswith(".mod"):
                    file_dupl = file.replace(".mod","_eff.mod")
                    shutil.copyfile(root+file,root+"eff\\"+file_dupl)
                #     with open(root+file) as f:
                #         lines_list = f.readlines()
                #         mod_files.append(lines_list)
                #         mod_names.append(file) #root+file   
    return mod_files, mod_names 

#read all .mod files in a directory and create the effective ones in parallel
def mod_to_eff(mech_folder, restrictions = None):
    mod_files = []
    mod_names = []
    for root, dirs, files in os.walk(mech_folder):
        for file in files:
            if restrictions:
                if file.endswith(".mod") and not file.endswith("_eff.mod") and file.replace(".mod","") in restrictions:
                    file_dupl = file.replace(".mod","_eff.mod")
                    with open (root+file) as f, open(root+"eff\\"+file_dupl,'w') as dupl:
                        for line in f:
                            if '{' in line:
                                dupl.write(line)
            elif not file.endswith("_eff.mod"):
                if file.endswith(".mod"):
                    file_dupl = file.replace(".mod","_eff.mod")
                    with open (root+file) as f, open(root+"eff\\"+file_dupl,'w') as dupl:
                        for line in f:
                            if '{' in line:
                                dupl.write(line)
    return mod_files, mod_names 

#read all all gbars from biophysics.hoc and put them into a dictionary
def read_gbars(cell_folder,d_2_s):

    #d_2_s = 0.001 #distance to soma is taken here as a constant but should be adapted to each section

    gbar_dict = {}
    lines_list = []
    for root, dirs, files in os.walk(cell_folder):
        for file in files:
            if file.endswith("biophysics.hoc"):
                with open(root+file) as f:
                    lines_list = f.readlines()
    for line in lines_list:
        search = re.search('g.*bar',line)
        if search:
            gbar = search.group().lower()
            location = re.search("\$o1\.[a-zA-Z]*",line).group().replace('$o1.','').lower()
            expression = re.search("\"(?![a-zA-Z,]).*\"\)",line).group().replace('%g','d_2_s')[1:-2].replace('exp','np.exp').replace('^','**').lower()
            if location in gbar_dict: 
                gbar_dict[location][gbar] = eval(expression)*1e4  #(S/cm2 -> S/m2)
            else:
                gbar_dict[location] = {gbar : eval(expression)*1e4} #(S/cm2 -> S/m2)
    return gbar_dict 

#filters in the list of the model file texts the lines that contain a gating parameter
def filter_mod(mod_files,mod_names):
    alpha_x = ".*[Aa]lpha.*=.*v.*|[Aa]lpha.*.*=.*v.*" # '.*' can be anything, '|': OR
    beta_x = ".*[Bb]eta.*=.*v.*|[Bb]eta.*.*=.*v.*" # '.*' can be anything, '|': OR
    tau_x = ".*[Tt]au.*=.*v.*|[Tt]au.*.*=.*v.*" # '.*' can be anything, '|': OR
    x_inf = ".*[Ii]nf.*=.*v.*|[Ii]nf.*.*=.*v.*" # '.*' can be anything, '|': OR

    l_alphas, l_betas, l_taus, l_infs = [], [], [], []
    hits = []
    for i,e in enumerate(mod_files):
        for j,f in enumerate(e):
            #print(test,f)
            l_alphas.append(([i,j],f,mod_names[i])) if  re.search(alpha_x,f) is not None else None
            l_betas.append(([i,j],f,mod_names[i])) if  re.search(beta_x,f) is not None else None
            l_taus.append(([i,j],f,mod_names[i])) if  re.search(tau_x,f) is not None else None
            l_infs.append(([i,j],f,mod_names[i])) if  re.search(x_inf,f) is not None else None
            hits.append(i) if  (re.search(alpha_x,f) or re.search(beta_x,f) or re.search(tau_x,f) or re.search(x_inf,f)) else None
    return l_alphas, l_betas, l_taus, l_infs, hits

#determine all the states from the 4 lists -> states names and descriptions
def states_from_lists(l_alphas, l_betas, l_taus, l_infs):

    alpha_pattern = "[a-zA-Z]_[Aa]lpha|[Aa]lpha_[a-zA-Z]|[a-zA-Z][Aa]lpha|[Aa]lpha[a-zA-Z]"
    beta_pattern = "[a-zA-Z]_[Bb]eta|[Bb]eta_[a-zA-Z]|[a-zA-Z][Bb]eta|[Bb]eta[a-zA-Z]"
    tau_pattern = "[a-zA-Z]_[Tt]au|[Tt]au_[a-zA-Z]|[a-zA-Z][Tt]au|[Tt]au[a-zA-Z]"
    inf_pattern = "[a-zA-Z]_[Ii]nf|[Ii]nf_[a-zA-Z]|[a-zA-Z][Ii]nf|[Ii]nf[a-zA-Z]"

    d_states = {}
    for e in l_alphas:
        match = re.search(alpha_pattern,e[1])
        if re.search("^h|h$",match.group()):
            key = "h_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" inactivation gate"
        elif re.search("^m|m$",match.group()):
            key = "m_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" activation gate"       
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        d_states[key] = value

    for e in l_betas:
        match = re.search(beta_pattern,e[1])
        if re.search("^h|h$",match.group()):
            key = "h_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" inactivation gate"
        elif re.search("^m|m$",match.group()):
            key = "m_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" activation gate"       
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        d_states[key] = value

    for e in l_taus:
        match = re.search(tau_pattern,e[1])
        if re.search("^h|h$",match.group()):
            key = "h_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" inactivation gate"
        elif re.search("^m|m$",match.group()):
            key = "m_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" activation gate"       
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        d_states[key] = value

    for e in l_infs:
        match = re.search(inf_pattern,e[1])
        if re.search("^h|h$",match.group()):
            key = "h_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" inactivation gate"
        elif re.search("^m|m$",match.group()):
            key = "m_"+e[2].replace('.mod','')
            value = e[2].replace('.mod','')+" activation gate"       
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        d_states[key] = value

    return d_states    

# extract both the LHS and RHS from the equations in the list
def formulas_from_lists(l_alphas, l_betas, l_taus, l_infs):

    alpha_pattern = "[a-zA-Z]_[Aa]lpha|[Aa]lpha_[a-zA-Z]|[a-zA-Z][Aa]lpha|[Aa]lpha[a-zA-Z]"
    beta_pattern = "[a-zA-Z]_[Bb]eta|[Bb]eta_[a-zA-Z]|[a-zA-Z][Bb]eta|[Bb]eta[a-zA-Z]"
    tau_pattern = "[a-zA-Z]_[Tt]au|[Tt]au_[a-zA-Z]|[a-zA-Z][Tt]au|[Tt]au[a-zA-Z]"
    inf_pattern = "[a-zA-Z]_[Ii]nf|[Ii]nf_[a-zA-Z]|[a-zA-Z][Ii]nf|[Ii]nf[a-zA-Z]"
    math_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z][ 0-9\.\+\-\*/\(\)a-zA-Z]*[0-9\.\+\-\*/\(\)a-zA-Z]" #removal of \t, \n and spaces around the formula

    d_alphas, d_betas, d_taus, d_infs = {}, {}, {}, {}
    for e in l_alphas:
        LHS,RHS = e[1].split("=")
        match = re.search(alpha_pattern,LHS)
        if re.search("^h|h$",match.group()):
            key = "h_"+'alpha_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        elif re.search("^m|m$",match.group()):
            key = "m_"+'alpha_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        value = value.replace('exp','np.exp')
        d_alphas[key] = value

    for e in l_betas:
        LHS,RHS = e[1].split("=")
        match = re.search(beta_pattern,LHS)
        if re.search("^h|h$",match.group()):
            key = "h_"+'beta_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        elif re.search("^m|m$",match.group()):
            key = "m_"+'beta_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        value = value.replace('exp','np.exp')
        d_betas[key] = value

    for e in l_taus:
        LHS,RHS = e[1].split("=")
        match = re.search(tau_pattern,LHS)
        if re.search("^h|h$",match.group()):
            key = "h_"+'tau_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        elif re.search("^m|m$",match.group()):
            key = "m_"+'tau_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        value = value.replace('exp','np.exp')
        d_taus[key] = value

    for e in l_infs:
        LHS,RHS = e[1].split("=")
        match = re.search(inf_pattern,LHS)
        if re.search("^h|h$",match.group()):
            key = "h_"+'inf_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        elif re.search("^m|m$",match.group()):
            key = "m_"+'inf_'+e[2].replace('.mod','')
            value = re.search(math_pattern,RHS).group()
        else:
           print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC}')
        value = value.replace('exp','np.exp')
        d_infs[key] = value
    return d_alphas, d_betas, d_taus, d_infs

#steady states are calculated based on alphas, betas and infs
def steadystates_from_gating_old(alphas, betas,taus,infs,dstates):
    #print(dstates)
    steadystates = {}
    for e in dstates:
        pre_suf = e.split('_',1) #only split on the first underscore
        if pre_suf[0]+'_alpha_'+pre_suf[1] in alphas.keys():
            alpha = alphas[pre_suf[0]+'_alpha_'+pre_suf[1]]
            beta = betas[pre_suf[0]+'_beta_'+pre_suf[1]]
            steadystates[e] = lambda Vm : eval(alpha)/(eval(alpha)+eval(beta))
        else:
            inf = infs[pre_suf[0]+'_inf_'+pre_suf[1]]
            steadystates[e] = lambda Vm : eval(inf)
    return steadystates

#steady states are calculated based on the values of alpha, beta and inf
def steadystates_from_gating(states,gating_states_kinetics):
    steadystates = {}
    for e in states:
        activ, mech = e.split('_',1)
        mech_var = gating_states_kinetics[mech]
        if activ+'alpha' in mech_var.keys():
            alpha = mech_var[activ+'alpha']
            beta = mech_var[activ+'beta']
            steadystates[e] = alpha/(alpha+beta)
        else:
            inf = mech_var[activ+'inf']    
            steadystates[e] = inf
    return steadystates

#derivative states functions are calculated based on the gating states kinetics and put into dictionary
def derstates_from_gating(states,gating_states_kinetics,x_dict):
    derstates = {}
    for e in states:
        activ, mech = e.split('_',1)
        mech_var = gating_states_kinetics[mech]
        if activ+'alpha' in mech_var.keys():
            mech_var[activ+'alpha']
            derstates[e] = mech_var[activ+'alpha'] * (1-x_dict[e]) - mech_var[activ+'beta'] * x_dict[e]
        else:
            mech_var[activ+'tau']
            derstates[e] = (mech_var[activ+'inf']-x_dict[e])/mech_var[activ+'tau']
    return derstates

#parses an equation and puts the excecuted RHS into the variable of the LHS
def calc_eq(e,var_pattern,math_pattern,equation_pattern,variables_dict):
    LHS,RHS = e.split("=")
    variable = re.search(var_pattern,LHS).group()
    formula = re.search(math_pattern,RHS).group()
    equation = re.search(equation_pattern,e).group()
    variable = variable.replace('exp','np.exp').replace('^','**').lower()
    formula = formula.replace('exp','np.exp').replace('^','**').lower()
    equation = equation.replace('exp','np.exp').replace('^','**').lower()
    try:
        variables_dict[variable] = eval(formula,{**globals(),**variables_dict}) # eval evaluates the value of the formula (RHS)
        exec(equation,{**globals(),**variables_dict}) # exec executes the equation and puts the value into the LHS
    except: 
        print(f"{tc.bcolors.OKCYAN}LOG: \t didn't work to compute: {equation}{tc.bcolors.ENDC}")

# a recursive functions which handels executing if-statements encountered in a MODL file
def if_recursive(list_mod,if_stat,variables_dict,offset):
    var_pattern = "[a-zA-Z_][a-zA-Z0-9_]*"
    math_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z][ 0-9\.\+\-\*/\(\)a-zA-Z^_]*[0-9\.\+\-\*/\(\)a-zA-Z]" #removal of \t, \n and spaces around the formula
    equation_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z][ 0-9\.\+\-\*/\(\)a-zA-Z^=_]*[0-9\.\+\-\*/\(\)a-zA-Z]" #removal of \t, \n and spaces around the formula
    if_pattern = "\(.*\)"

    while not re.search("}",list_mod[offset]):
        if re.search('if',list_mod[offset]) and if_statement:
            if_statement = eval(re.search(if_pattern,list_mod).group().replace('exp','np.exp').replace('^','**').lower(),{**globals(),**variables_dict})
            if_recursive(list_mod[offset+1],if_statement,variables_dict,offset+1)
        elif re.search('else',list_mod[i]) and if_statement:
            return variables_dict, offset+1
        elif re.search('=',e) and if_statement:
            calc_eq(list_mod[offset],var_pattern,math_pattern,equation_pattern,variables_dict)
            offset += 1   
    return variables_dict, offset+1    

#extract PROCEDURE block and execute all equations to retrieve alpha, beta, tau and inf from m and h
def gating_from_PROCEDURES(list_mod,mod_name,Vm,start_executing = 0):
    var_pattern = "[a-zA-Z_][a-zA-Z0-9_]*"
    math_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z][ 0-9\.\+\-\*/\(\)a-zA-Z^_]*[0-9\.\+\-\*/\(\)a-zA-Z]" #removal of \t, \n and spaces around the formula
    equation_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z][ 0-9\.\+\-\*/\(\)a-zA-Z^=_]*[0-9\.\+\-\*/\(\)a-zA-Z]" #removal of \t, \n and spaces around the formula
    if_pattern = "\(.*\)"
    celsius = tc.T_C #temperature in degrees Celsius

    variables_dict = {'v': Vm, 'celsius' : celsius} #uncomment this if v is a known parameter
    i = 0
    while i < len (list_mod):
        e = list_mod[i]
    #for i,e in enumerate(list_mod):
        #stop when PROCEDURE BLOCK is finished
        if re.search('^}',e):
            start_executing = 0
        if start_executing == 1:
            if re.search('if',e):
                if_statement = eval(re.search(if_pattern,e).group().replace('exp','np.exp').replace('^','**').lower(),{**globals(),**variables_dict})
                variables_dict,offset = if_recursive(list_mod[i+1:],if_statement,variables_dict,1)
                i += offset
                continue
            #elif?
            if re.search('=',e):
                calc_eq(e,var_pattern,math_pattern,equation_pattern,variables_dict)
        #start when entering PROCEDURE block on the next iteration
        if re.search('PROCEDURE rate',e):
            start_executing = 1
        i += 1

    #maybe good idea to add a suffix to each variable -> NO
    #variables_dict[variable+'_'+mod_name.replace(".mod","")]
    for x in ['m','h']:
        if not x+'alpha' in variables_dict.keys() and x+'tau' in variables_dict.keys() and x+'inf' in variables_dict.keys():
            variables_dict[x+'alpha'] =  variables_dict[x+'inf'] / variables_dict[x+'tau']
        if not x+'beta' in variables_dict.keys() and x+'tau' in variables_dict.keys() and x+'inf' in variables_dict.keys():
            variables_dict[x+'beta'] =  variables_dict[x+'inf'] / variables_dict[x+'tau']
    return variables_dict   

#extract PROCEDURE block and execute all equations to retrieve alpha, beta, tau and inf from m and h
def currents_from_BREAKPOINT(list_mod,mod_name,Vm,x_dict,g_dict,location,start_executing = 0):
    var_pattern = "[a-zA-Z_][a-zA-Z0-9_]*"
    math_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z][ 0-9\.\+\-\*/\(\)a-zA-Z^_]*[0-9\.\+\-\*/\(\)a-zA-Z]" #removal of \t, \n and spaces around the formula
    equation_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z][ 0-9\.\+\-\*/\(\)a-zA-Z^=_]*[0-9\.\+\-\*/\(\)a-zA-Z]" #removal of \t, \n and spaces around the formula
    if_pattern = "\(.*\)"
    celsius = tc.T_C #temperature in degrees Celsius
    #gbar = eval('g'+mod_name.replace('.mod','')+'bar')
    ek, ena, eleak, ehcn = tc.EK, tc.ENa, tc.ELeak, tc.ehcn

    try:
        m = x_dict['m_'+mod_name.replace('.mod','')]
    except:
        m = 0
        print(f'{tc.bcolors.OKCYAN}LOG: \t {x_dict} does not contain an m for: {mod_name}{tc.bcolors.ENDC}')
    try:
        h = x_dict['h_'+mod_name.replace('.mod','')]
    except:
        h = 0
        print(f'{tc.bcolors.OKCYAN}LOG: \t {x_dict} does not contain an h for {mod_name}{tc.bcolors.ENDC}')

    variables_dict = {'v': Vm, 'celsius' : celsius, 'm' : m,\
                      'h' : h, 'ek' : ek, 'ena' : ena,\
                        'eleak' : eleak, 'ehcn' : ehcn} | g_dict[location]
    variables_dict_backup = copy.copy(variables_dict)
    i = 0
    while i < len (list_mod):
        e = list_mod[i]
    #for i,e in enumerate(list_mod):
        #stop when PROCEDURE BLOCK is finished
        if re.search('^}',e):
            start_executing = 0
        if start_executing == 1:
            if re.search('if',e):
                if_statement = eval(re.search(if_pattern,e).group().replace('exp','np.exp').replace('^','**').lower(),{**globals(),**variables_dict})
                variables_dict,offset = if_recursive(list_mod[i+1:],if_statement,variables_dict,1)
                i += offset
                continue
            #elif?
            if re.search('=',e):
                calc_eq(e,var_pattern,math_pattern,equation_pattern,variables_dict)
        #start when entering PROCEDURE block on the next iteration
        if re.search('BREAKPOINT',e):
            start_executing = 1
        i += 1

    #maybe good idea to add a suffix to each variable
    #variables_dict[variable+'_'+mod_name.replace(".mod","")]

    # reduced dictionary with only the newly calculated values
    variables_reduced = copy.copy(variables_dict)
    for e in variables_dict_backup:
        del variables_reduced[e]
    return variables_dict, variables_reduced

#reduce the number of underscores in a name/string to 1
def rm_us(e): 
    while e.count('_') > 0:
        suf,pre = e[::-1].split('_',1) # reverse the string and remove the first underscore
        e = pre[::-1]+suf[::-1] #re-reverse the string in order that the last suffix is removed
    return e