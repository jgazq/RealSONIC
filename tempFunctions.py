""""contains all the utility/auxiliary functions for both NEURONic as Pythonic"""


"""-----------------------------------------------------------------------------------IMPORTS-----------------------------------------------------------------------------------"""
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
import zipfile

import tempConstants as tc
import prev.Interp3Dfield as tt


"""-----------------------------------------------------------------------------------UTILS-----------------------------------------------------------------------------------"""
def read_pickle(path, filename,prints=False):
    """read the LUT pickle files metadata"""

    pkl = pd.read_pickle(path+filename)
    if not prints:
        return pkl
    print('shape of V in tables',pkl['tables']['V'].shape) #dimensions of refs
    print('refs(inputs): ',pkl['refs'].keys()) #which variables are sweeped over
    #print('refs: ',pkl['refs'].values()) #the arrays of each variable(length of each variable is given in the dimensions above)
    print('tables(outputs): ',pkl['tables'].keys()) #which gating state parameters are stored in the LUT
    if 'tcomp' in pkl['tables'].keys(): #test files don't have a computation time value -> NOT TRUE!!! I ?accidentaly? removed the tcomp line in run_lookups_MC.py
        print('total simulation time: ',np.sum(pkl['tables']['tcomp'])/3600/24,'days')

    return pkl


def rm_us(name): 
    """reduce the number of underscores in a name/string to none"""
    name_ = []
    if type(name)==list:
        print('removing the underscores in a list')
        for e in name:
            while e.count('_') > 0:
                suf,pre = e[::-1].split('_',1) # reverse the string and remove the first underscore
                e = pre[::-1]+suf[::-1] #re-reverse the string in order that the last suffix is removed
            name_.append(e)
    else:
        while name.count('_') > 0:
            suf,pre = name[::-1].split('_',1) # reverse the string and remove the first underscore
            name = pre[::-1]+suf[::-1] #re-reverse the string in order that the last suffix is removed      
    return e if len(name_) > 1 else name


def unzip():
    """unzips/extracts all files in a certain folder -> here the folder of the cells of BBP is unzipped"""

    folder = "OneDrive - UGent/PhD/Untouched Code/BBP/zipped (1035 cells)"

    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            #print(file)
        # Path to the zip file
            zip_file_path = file#'OneDrive - UGent/PhD/Untouched Code/BBP/zipped (1035 cells)/L1_DAC_bNAC219_1.zip'

            # Directory where you want to extract the contents
            extract_dir = 'OneDrive - UGent/PhD/Untouched Code/BBP/unzipped'

            # Create the extract directory if it doesn't exist
            if not os.path.exists(extract_dir):
                os.makedirs(extract_dir)

            # Open the zip file
            with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
                # Extract all contents into the extract directory
                zip_ref.extractall(extract_dir)

            #print("Extraction complete!")

def allmech_BBP():
    """searches for all mechanism files and returns the intersection of these mechanisms in a list"""

    folder = "OneDrive - UGent/PhD/Untouched Code/BBP/unzipped"
    lijst = []
    locations = []

    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            filename = os.path.basename(file)
            #print(filename)
            if filename.endswith(".mod") and filename not in lijst:
                #print(lijst)
                lijst.append(filename)
                locations.append(file.split('unzipped')[-1])
    print(lijst,locations)
  

"""-----------------------------------------------------------------------------------INSIDE NEURON (MECHANISMS AND SECTIONS)-----------------------------------------------------------------------------------"""
def gbar_mechs():
    """get the conductance bar values for every mechanism, the value is taken from the first encounter"""

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


def loaded_mechs():
    """returns a list with all mechanisms that are loaded into hoc"""

    lijst = []
    mt = h.MechanismType(0)
    mname  = h.ref('')
    for i in range(mt.count()):
        mt.select(i)
        mt.selected(mname)
        lijst.append(mname[0])

    return lijst

                
def cell_mechs():
    """print all mechanisms of a cell """

    mech_actual = []
    mech_poss = loaded_mechs() #mech_poss contains all possible mechanisms that this cell can have
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

                  
def sec_mechs():
    """print every section with its corresponding mechanisms (every section contains segments
       that have the same mechanisms but are just spatially extended)"""

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


def dist_sec_soma():
    """prints the distance to the soma for all segments"""

    i = 0
    for sec in h.allsec():
        for seg in sec:
            if "Scale" not in str(sec) and "Elec" not in str(sec):
                x_SOM,y_SOM,z_SOM = tt.findSOMA()
                dist_2_soma = np.sqrt((seg.x_xtra-x_SOM)**2 + (seg.y_xtra-y_SOM)**2 + (seg.z_xtra-z_SOM)**2)
                print(seg,dist_2_soma)


def coord_dict():
    """returns a dict with all segments containing their 3D position"""

    i = 0
    dict = {}
    for sec in h.allsec():
        for seg in sec:
            if "Scale" not in str(sec) and "Elec" not in str(sec):
                dict[str(seg)] = (seg.x_xtra, seg.y_xtra, seg.z_xtra)

    return dict


def coord_dict_type():
    """returns a dict with each type of segment containing all the segments and their locations"""

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


"""-----------------------------------------------------------------------------------WRITE REALNEURON (AND MODL)-----------------------------------------------------------------------------------"""
def read_mod(mech_folder, restrictions = None):
    """read all .mod files in a directory and put them into a 2D-list
    restrictions: list that contains the mechanisms that need to be read"""

    mod_files = []
    mod_names = []
    for root, dirs, files in os.walk(mech_folder):
        for file in files:
            if restrictions: #list with restricted mechanisms is given
                if file.endswith(".mod") and not 'eff' in root+file and file.replace(".mod","") in restrictions:
                    with open(os.path.join(root,file)) as f:
                        lines_list = f.readlines()
                        mod_files.append(lines_list)
                        mod_names.append(file.replace('.mod','')) #root+file   
            elif not 'eff' in root+file: #no restriction list given
                if file.endswith(".mod"):
                    with open(os.path.join(root,file)) as f:
                        lines_list = f.readlines()
                        mod_files.append(lines_list)
                        mod_names.append(file.replace('.mod','')) #root+file  

    return mod_files, mod_names 


def mod_duplicate(mech_folder, restrictions = None):
    """read all .mod files in a directory and duplicate them in a seperate folder: eff"""
    mod_files = []
    mod_names = []
    for root, dirs, files in os.walk(mech_folder):
        for file in files:
            if restrictions:
                if file.endswith(".mod") and not file.endswith("_eff.mod") and (file.replace(".mod","") in restrictions or rm_us(file.replace(".mod","")) in restrictions):
                    file_dupl = file.replace(".mod","_eff.mod") #create a new file with _eff at the end
                    shutil.copyfile(os.path.join(root,file),os.path.join(root,"eff",file_dupl)) #and copy the original file into the newly created file
                else:
                    file_dupl = file #create a new file with _eff at the end
                    shutil.copyfile(os.path.join(root,file),os.path.join(root,"eff",file_dupl)) #and copy the original file into the newly created file                    
            elif not file.endswith("_eff.mod"):
                if file.endswith(".mod"):
                    file_dupl = file.replace(".mod","_eff.mod")
                    shutil.copyfile(os.path.join(root,file),os.path.join(root,"eff",file_dupl))

    return mod_files, mod_names 


def read_gbars(cell_folder,d_2_s):
    """read all all gbars from biophysics.hoc and put them into a dictionary"""

    #d_2_s = 0.001 #distance to soma is taken here as a constant but should be adapted to each section

    gbar_dict = {}
    lines_list = []
    for root, dirs, files in os.walk(cell_folder):
        for file in files:
            if file.endswith("biophysics.hoc"):
                with open(os.path.join(root,file)) as f:
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


def filter_mod(mod_files,mod_names):
    """filters in the list of the model file texts the lines that contain a gating parameter"""

    l_alphas, l_betas, l_taus, l_infs = [], [], [], []
    hits = []
    for i,e in enumerate(mod_files):
        for j,f in enumerate(e):
            #print(test,f)
            l_alphas.append(([i,j],f,mod_names[i])) if  re.search(tc.alpha_x_pattern,f) is not None else None
            l_betas.append(([i,j],f,mod_names[i])) if  re.search(tc.beta_x_pattern,f) is not None else None
            l_taus.append(([i,j],f,mod_names[i])) if  re.search(tc.tau_x_pattern,f) is not None else None
            l_infs.append(([i,j],f,mod_names[i])) if  re.search(tc.x_inf_pattern,f) is not None else None
            hits.append(i) if  (re.search(tc.alpha_x_pattern,f) or re.search(tc.beta_x_pattern,f) or re.search(tc.tau_x_pattern,f) or re.search(tc.x_inf_pattern,f)) else None

    return l_alphas, l_betas, l_taus, l_infs, hits

def state_from_line(line,pattern):
    """look if a gating state parameter is defined in the current/given line"""

    match = re.search(pattern,line[1])
    if re.search("^h|h$",match.group()):
        key = "h_"+line[2].replace('.mod','')
        value = line[2].replace('.mod','')+" inactivation gate"
    elif re.search("^m|m$",match.group()):
        key = "m_"+line[2].replace('.mod','')
        value = line[2].replace('.mod','')+" activation gate"       
    else:
        print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h in {line[2]}{tc.bcolors.ENDC}')
        return
    
    return key,value


def states_from_lists(l_alphas, l_betas, l_taus, l_infs):
    """determine all the states from the 4 lists -> states names and descriptions"""

    d_states = {}

    for e in l_alphas:
        if state_from_line(e,tc.alpha_pattern):
            key,value = state_from_line(e,tc.alpha_pattern)
            d_states[key] = value

    for e in l_betas:
        if state_from_line(e,tc.beta_pattern):
            key,value = state_from_line(e,tc.beta_pattern)
            d_states[key] = value

    for e in l_taus:
        if state_from_line(e,tc.tau_pattern):
            key,value = state_from_line(e,tc.tau_pattern)
            d_states[key] = value

    for e in l_infs:
        if state_from_line(e,tc.inf_pattern):
            key,value = state_from_line(e,tc.inf_pattern)
            d_states[key] = value

    return d_states    


def formula_from_line(line, pattern, kinetic):
    "split an equations into LHS and RHS"
    LHS,RHS = line[1].split("=")
    match = re.search(pattern,LHS)
    if re.search("^h|h$",match.group()):
        key = "h_"+kinetic+'_'+line[2].replace('.mod','')
        value = re.search(tc.math_pattern,RHS).group()
    elif re.search("^m|m$",match.group()):
        key = "m_"+kinetic+'_'+line[2].replace('.mod','')
        value = re.search(tc.math_pattern,RHS).group()
    else:
        print(f'{tc.bcolors.OKCYAN}{match} is not a recognized gating parameter m or h{tc.bcolors.ENDC} in {line[2]}')
        return
    value = value.replace('exp','np.exp') 

    return key,value   


def formulas_from_lists(l_alphas, l_betas, l_taus, l_infs):
    """extract both the LHS and RHS from the equations in the list"""

    d_alphas, d_betas, d_taus, d_infs = {}, {}, {}, {}
    for e in l_alphas:
        if formula_from_line(e,tc.alpha_pattern,"alpha"):
            key, value = formula_from_line(e,tc.alpha_pattern,"alpha")
            d_alphas[key] = value

    for e in l_betas:
        if formula_from_line(e,tc.beta_pattern,"beta"):
            key, value = formula_from_line(e,tc.beta_pattern,"beta")
            d_betas[key] = value

    for e in l_taus:
        if formula_from_line(e,tc.tau_pattern,"tau"):
            key, value = formula_from_line(e,tc.tau_pattern,"tau")
            d_taus[key] = value

    for e in l_infs:
        if formula_from_line(e,tc.inf_pattern,"inf"):
            key, value = formula_from_line(e,tc.inf_pattern,"inf")
            d_taus[key] = value

    return d_alphas, d_betas, d_taus, d_infs


def steadystates_from_gating_old(alphas, betas,taus,infs,dstates):
    """steady states are calculated based on alphas, betas and infs"""

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


def steadystates_from_gating(states,gating_states_kinetics):
    """steady states are calculated based on the values of alpha, beta and inf"""

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


def derstates_from_gating(states,gating_states_kinetics,x_dict):
    """derivative states functions are calculated based on the gating states kinetics and put into dictionary"""

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


def eq_hoc2pyt(equation):
    "translate an equation or formula from hoc syntax to python syntax"

    equation = equation.replace('exp','np.exp').replace('^','**')
    equation = equation.strip() #remove spaces at beginning and end
    equation = equation.lower() #put equations always in lowercase #POTENTIAL_RISK

    return equation


def calc_eq(e,variables_dict,mod_name=None):
    """parses an equation and puts the excecuted RHS into the variable of the LHS"""

    LHS,RHS = e.split("=")
    variable = re.search(tc.var_pattern,LHS).group()
    formula = re.search(tc.math_pattern,RHS).group()
    equation = re.search(tc.equation_pattern,e).group()
    variable = eq_hoc2pyt(variable)
    formula = eq_hoc2pyt(formula)
    equation = eq_hoc2pyt(equation)
    # try:
    # if mod_name == 'Ca':
    #     print(f'equation: {variable} : {formula} : {equation}')
    if True:
        variables_dict[variable] = eval(formula,{**globals(),**variables_dict}) # eval evaluates the value of the formula (RHS)
        exec(equation,{**globals(),**variables_dict}) # exec executes the equation and puts the value into the LHS
    # except: 
    #     print(f"{tc.bcolors.OKCYAN}LOG: \t didn't work to compute: {equation}{tc.bcolors.ENDC}")
    #     print(equation,variables_dict)
    #     quit()


def if_recursive(list_mod,if_statement,variables_dict,offset,mod_name=None):
    """a recursive functions which handels executing if-statements encountered in a MODL file"""
    orig_offset = offset
    # if mod_name == 'Ca':
    #     print('if recursive called')
    #     print("list_mod:", list_mod)
    while not re.search("[^a-zA-Z]}[^a-zA-Z]",list_mod[offset-orig_offset]):
        # if mod_name == 'Ca':
        #     print(f'list_mod[offset-orig_offset]: {list_mod[offset-orig_offset]}')
        if re.search('if',list_mod[offset-orig_offset]) and if_statement:
            # if mod_name == 'Ca':
            #     print('here1')
            if_statement = eval(re.search(tc.if_pattern,list_mod).group().replace('exp','np.exp').replace('^','**').lower(),{**globals(),**variables_dict})
            if_recursive(list_mod[offset+1-orig_offset:],if_statement,variables_dict,offset+1,mod_name)
        elif re.search('else',list_mod[offset-orig_offset]): #and if_statement:
            # if mod_name == 'Ca':
            #     print('here2')
            return if_recursive(list_mod[offset+1-orig_offset:],not if_statement,variables_dict,offset+1,mod_name) #return variables_dict, offset+1
        elif re.search('=',list_mod[offset-orig_offset]) and if_statement:
            # if mod_name == 'Ca':
            #     print('here3')
            #     print(list_mod[offset-orig_offset])
            #     print(variables_dict)
            calc_eq(list_mod[offset-orig_offset],variables_dict,mod_name=mod_name)
            offset += 1 
        else:
            # if mod_name == 'Ca':
            #     print('here4')
            offset +=1  

    return variables_dict, offset+1    


def gating_from_PROCEDURES(list_mod,mod_name,Vm): #,start_executing = 0):
    """extract PROCEDURE block and execute all equations to retrieve alpha, beta, tau and inf from m and h"""

    #print(mod_name)
    celsius = tc.T_C #temperature in degrees Celsius
    proc_executing, param_executing = 0,0

    variables_dict = {'v': Vm, 'celsius' : celsius} #uncomment this if v is a known parameter
    # if mod_name == 'Ca':
    #     print(f"voltage = {variables_dict['v']}")
    i = 0
    while i < len (list_mod):
        e = list_mod[i]
        # if i > 55 and mod_name == 'Ca':
        #     print("line: ",e)
    #for i,e in enumerate(list_mod):
        #stop when PROCEDURE BLOCK is finished
        if re.search('^}',e):
            proc_executing = 0
            param_executing = 0
        if proc_executing == 1:
            if re.search('if',e):
                # if mod_name == 'Ca':
                #     variables_dict["v"] = -27
                #     print(f'mod_name: {mod_name}, v = {variables_dict["v"]} , line: {e}')
                if_statement = eval(re.search(tc.if_pattern,e).group().replace('exp','np.exp').replace('^','**').lower(),{**globals(),**variables_dict})
                # if mod_name == 'Ca':
                #     print('if_statement: ',if_statement)
                variables_dict,offset = if_recursive(list_mod[i+1:],if_statement,variables_dict,0,mod_name)
                # if mod_name == 'Ca':
                #     print(f"variables_dict : {variables_dict}, offset: {offset}")                    
                i += offset
                continue
            #elif?
            if re.search('=',e):
                calc_eq(e,variables_dict,mod_name)
        elif param_executing == 1:
            if re.search('=',e):
                if re.findall("\(.*\)",e):
                    calc_eq(e.replace(re.findall("\(.*\)",e)[0],""),variables_dict) #when evaluating the equations in PARAM block, the units are removed that are between brackets  
                else:
                    calc_eq(e,variables_dict) #don't replace anything if units are not given
        #start when entering PROCEDURE block on the next iteration
        if re.search('PROCEDURE',e):
            proc_executing = 1
        #start when entering PARAMETER block on the next iteration
        if re.search('PARAMETER',e):
            param_executing = 1
        i += 1

    #maybe good idea to add a suffix to each variable -> NO
    #variables_dict[variable+'_'+mod_name.replace(".mod","")]
    for x in ['m','h']:
        if not x+'alpha' in variables_dict.keys() and x+'tau' in variables_dict.keys() and x+'inf' in variables_dict.keys():
            variables_dict[x+'alpha'] = variables_dict[x+'inf'] / variables_dict[x+'tau']
        if not x+'beta' in variables_dict.keys() and x+'tau' in variables_dict.keys() and x+'inf' in variables_dict.keys():
            variables_dict[x+'beta'] = (1 - variables_dict[x+'inf']) / variables_dict[x+'tau']
    # if mod_name == 'Ca':
    #     print("final variables_dict: ",variables_dict)
    #     quit()
    return variables_dict   


def currents_from_BREAKPOINT(list_mod,mod_name,Vm,x_dict,g_dict,location,start_executing = 0):
    """extract PROCEDURE block and execute all equations to retrieve alpha, beta, tau and inf from m and h"""

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
                if_statement = eval(re.search(tc.if_pattern,e).group().replace('exp','np.exp').replace('^','**').lower(),{**globals(),**variables_dict})
                variables_dict,offset = if_recursive(list_mod[i+1:],if_statement,variables_dict,1)
                i += offset
                continue
            #elif?
            if re.search('=',e):
                calc_eq(e,variables_dict)
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


"""-----------------------------------------------------------------------------------WRITE MODL-----------------------------------------------------------------------------------"""
def one_to_multiline(root,file):
    """split a BLOCK defined a .modl file on one line over multiple lines {
    like this
    }"""
    changed = 0
    with open(os.path.join(root,file), 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if re.search(tc.onelineBLOCK_pattern,line):
            part1, part2 = line.split("{")[0],line.split("{")[1].split("}")[0]
            lines[i] = f"{part1}{{\n\t{part2}\n}}"
            (line.split("{")[0],"\n",line.split("{")[1].split("}")[0])
            changed = 1
            print(line,file)
            #break
    if changed:
        with open(os.path.join(root,file), 'w') as newf: #.replace(".mod","_test.mod")
            newf.writelines(lines)

        return file #.replace(".mod","_test.mod")
    else:

        return file
    

def model_to_BLOCKS(flist):
    """put the different BLOCK lines of a mechanism .modl file in a dictionary with each seperate block"""

    BLOCKdict = {}
    BLOCK_ON = 0
    BLOCK_list = []
    for line in flist:
        if re.search(tc.block_init_pattern,line): #do determine if we are in a specific block or not
            block = re.search(tc.block_pattern,(re.search(tc.block_init_pattern,line).group(0))).group(0) #first look if we are in a BLOCK initiation line and then extract the actual block
            BLOCK_ON = 1
        if BLOCK_ON == 1:
            BLOCK_list.append(line)
        if line.startswith('}'):
            BLOCK_ON = 0
            BLOCKdict[block] = BLOCK_list
            BLOCK_list = []
    return BLOCKdict


def str_in_element_in_list(str,list):
    """looks if a string is inside an element of a list(to search for certain keywords)"""
    return [e for e in list if str in e]


def str1_before_str2(str1,str2,flist):
    """str1 NEEDS to be before str2 so interchange strings if str2 is before str1"""
    BLOCKdict = model_to_BLOCKS(flist)
    for e,f in BLOCKdict.items():
        list1 = str_in_element_in_list(str1,f) #local
        list2 = str_in_element_in_list(str2,f) #update
        if list1 and list2:
            loc1 = f.index(list1[0])
            loc2 = f.index(list2[0])
            if loc1 > loc2: #interchange them if str2 is before str1

                return f[loc1],f[loc2]

        # if any([str1 in g for g in f]) and any([str2 in g for g in f]):
        #     str1_index = f.index(str1)
        #     str2_index = f.index(str1)


            

"""-----------------------------------------------------------------------------------VARIA-----------------------------------------------------------------------------------"""
def plt_transdistr(psource,grid_type):
    """plot the distribution points of a transducer source"""

    nsources = psource.get_default_nsources()
    x, y = psource.getXYSources(m=nsources, d=grid_type)
    plt.title(f'{grid_type} distribution')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.plot(x * ms.M_TO_MM, y * ms.M_TO_MM, 'o', markersize = 4)
    plt.tight_layout()
    plt.show()


def plt_LUTeff(variable,table,Q_refs,A_refs,plot=False,reduced_range=False):
    """"plot the given effective variable in function of the charge density for different pressure amplitudes
    this function assumes that the LUT is calculated for only 1 sonophore radius, 1 frequency and 1 sonophore coverage fraction"""
    num_shades = len(Q_refs)

    cmap = plt.cm.get_cmap('hsv')
    # Generate colors across the colormap
    colors = [cmap(i / num_shades) for i in range(num_shades)]
    #colors = plt.cm.Reds(np.linspace(0, 1, len(pkldict['refs']['Q'])))
    plt.clf()
    no_reduc = ['tcomp','V']
    for i,e in enumerate(table[::5]):
        plt.plot(Q_refs[::5]*1e5,e[::5],label=f"{A_refs[::5][i]*1e-3:.1e} kPa".replace('+0','').replace('-0',''),color = colors[::5][i])
        plt.legend()
        plt.xlabel('$\mathrm{Q [nC/cm^2]}$')
        ylab = 'V [mV]' if variable == 'V' else 'C $\mathrm{[\dfrac{uF}{cm^2}]}$' if variable == 'C' else variable+ ' $\mathrm{[\dfrac{1}{ms}]}$' if ('alpha' in variable or 'beta' in variable) else variable
        plt.ylabel(f'{ylab}')
        if reduced_range and (variable not in no_reduc):
            change_lower = plt.ylim()[0] < -10
            change_upper = plt.ylim()[1] > 10
            e_sorted = np.sort(e[::5].reshape(-1))
            e_limited = [e for e in e_sorted if e < 10 and e > -10]
            if not e_limited: #if no value falls under the limited range, use the original array instead of no values
                e_limited = e
            # print(variable); print('before'); print(plt.ylim())
            plt.ylim(bottom = min(np.min(e_limited),-0.5)) if change_lower else None
            plt.ylim(top = max(np.max(e_limited),0.5)) if change_upper else None
            # print('after'); print(plt.ylim()); print('\n\n')
        if plot:
            plt.show()
        


def save_gatingplots(pkldict,foldername,reduced_range=True):
    """"save all plots for the different LUTs, containing the effective variables of the gating kinetics
    (this function assumes that the LUT is calculated for only 1 sonophore radius, 1 frequency and 1 sonophore coverage fraction)"""

    Qrange,Arange = pkldict['refs']['Q'], pkldict['refs']['A']
    plt.figure(figsize=(8, 6)) #change figsize so all plots are shown properly and fit in the box
    #plot all calculated effective variables
    for key in pkldict['tables']:
        table = pkldict['tables'][key]
        table = table[0][0]
        plt_LUTeff(key,table,Qrange,Arange,reduced_range=reduced_range)
        plt.savefig(f'figs/{foldername}/{key}.png')

    #plot the effective capacity
    table_V = pkldict['tables']['V'][0][0]
    Q_table = np.tile(Qrange,np.prod(table_V.shape)//len(Qrange)).reshape(table_V.shape)
    table_C = Q_table / table_V
    plt_LUTeff('C',table_C,Qrange,Arange)
    plt.savefig(f'figs/{foldername}/C.png')   
