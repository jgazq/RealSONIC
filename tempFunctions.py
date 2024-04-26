""""contains all the utility/auxiliary functions for both NEURONic as Pythonic"""


"""-----------------------------------------------------------------------------------IMPORTS-----------------------------------------------------------------------------------"""
import matplotlib.pyplot as plt
# import PySONIC as ps
# import MorphoSONIC as ms
import numpy as np
import pickle
#import pandas as pd
from neuron import h
import re
import os
import copy
import shutil
import zipfile
import datetime
import scipy.interpolate as interp

import tempConstants as tc
import prev.Interp3Dfield as tt


"""-----------------------------------------------------------------------------------UTILS-----------------------------------------------------------------------------------"""
def read_pickle(path, filename=None, prints=False):
    """ read the LUT pickle files metadata
        :path: directory where the pickle file is stored (can contain also the filename)
        :filename: optional to give the filename separate 
        :prints: prints various metadata of the loaded LUT"""

    fpath = path+filename if filename else path #option to give path and filename together or seperate
    #pkl = pd.read_pickle(path+filename)
    with open(fpath, 'rb') as fh:
        pkl = pickle.load(fh)
    if not prints:
        return pkl
    print('shape of V in tables',pkl['tables']['V'].shape) #dimensions of refs
    print('refs(inputs): ',pkl['refs'].keys()) #which variables are sweeped over
    print('refs values (inputs): ',pkl['refs'].values())
    #print('refs: ',pkl['refs'].values()) #the arrays of each variable(length of each variable is given in the dimensions above)
    print('tables(outputs): ',pkl['tables'].keys()) #which gating state parameters are stored in the LUT
    if 'tcomp' in pkl['tables'].keys(): #test files don't have a computation time value -> NOT TRUE!!! I ?accidentaly? removed the tcomp line in run_lookups_MC.py
        print('total simulation time: ',np.sum(pkl['tables']['tcomp'])/3600/24,'days')
    # for e,f in pkl['tables'].items(): #to print all the tables for all gating parameters
    #     print(f'{e}:\n\n{f}')
    return pkl


def load_pickle(dictio, path, filename = None):
    """ saves a dictionary into a pickle file with the provided path and name
        :dictio: dictionary containing the LUT
        :path: directory where the LUT needs to be stored (filename can also be included in path)
        :filename: explicit mentioning the name of the file where it needs to be stored"""

    fpath = path+filename if filename else path #option to give path and filename together or seperate
    with open(fpath, 'wb') as fh:
        pickle.dump(dictio, fh)
    print(f'loaded pickle file with shape: {dictio["tables"]["V"].shape} in: {fpath}')
    

def rm_us(name): 
    """ reduce the number of underscores in a name/string to none
        :name: string which needs underscore reduction"""

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

    folder = r"C:\Users\jgazquez\OneDrive - UGent\PhD\Untouched Code\BBP\unzipped"
    lijst = []
    locations = []
    mech_dict = {}
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            filename = os.path.basename(file)
            # print(filename)
            if filename.endswith(".mod") and filename in mech_dict.keys():
                mech_dict[filename].append(file.split('\\')[-3])
            elif filename.endswith(".mod"):
                mech_dict[filename] = [file.split('\\')[-3]]
            if filename.endswith(".mod") and filename not in lijst:
                #print(lijst)
                lijst.append(filename)
                locations.append(file.split('unzipped')[-1])
    #print(f"{'StochKv.mod'} found in the following cells:\n{mech_dict['StochKv.mod']}\n\n")
    print(lijst,locations)


def tcomp_run_time(tcomp,run_time):
    """ divides the tcomp value from PySONIC by the actual run_time to estimate the number of used cores
        :tcomp: serial computation time, format: H:M:S.MS s
        :run_time: actual, parallel computation time,format: x.y days"""
    if ',' in run_time:
        days, run_time = run_time.split(',')
        days = re.findall("\d*",days)[0] #remove unit
    else:
        days = 0
    hours,ms = run_time.split('.')
    ms = re.findall("\d*",ms)[0] #remove unit
    h,m,s = hours.split(':')
    run_time_days = float(days) + float(h)/24 + float(m)/24/60 + float(s)/24/60/60 + float(ms)/24/60/60*1e-3
    tcomp_days = float(re.findall("\d*\.\d*",tcomp)[0]) #remove unit
    return tcomp_days/run_time_days


def write_csv(numA=50, Qstep=10, maxovertones=1, Vm0=-75, Cm0=2e-2, max_lines=np.inf):
    """ writes all the parameter combinations into a csv file which can be used with the worker module for parallelisation
        :numA: number of amplitude values
        :Qstep: the difference between two charge values
        :maxovertones: the maximum number of overtones
        :Vm0: the lowest possible membrane resting voltage (worst case)
        :Cm0: the highest possible membrane resting capacitance (worst case)
        :max_lines: the maximum number of lines that can be written out to the .csv file"""

    nlines = 0
    DQ_LOOKUP = 1e-5
    numQ = int(((50.0 - np.round(Vm0 - 35.0)) * Cm0 * 1e-3) / DQ_LOOKUP + 1)
    with open("lookup_data_ext.csv",'w') as file:
        file.write("radius, frequency, amplitude, startcharge, endcharge, overtones\n")
        for overt in range(maxovertones+1):
            for a in np.array([16.0, 32.0, 64.0]):
                for f in np.array([20., 100., 500., 1e3, 2e3, 3e3, 4e3]):
                    for P_A in np.insert(np.logspace(np.log10(0.1), np.log10(600), num=numA), 0, 0.0):
                        for Qstart in range(0,numQ,Qstep):
                            file.write(f"{a}, {f}, {P_A}, {Qstart}, {Qstart+Qstep}, {overt}\n")
                            nlines += 1
                            if nlines == max_lines:
                                print('Maximum number of lines reached')
                                quit()


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
    dictio = {}
    for sec in h.allsec():
        for seg in sec:
            if "Scale" not in str(sec) and "Elec" not in str(sec):
                dictio[str(seg)] = (seg.x_xtra, seg.y_xtra, seg.z_xtra)

    return dictio


def coord_dict_type():
    """returns a dict with each type of segment containing all the segments and their locations"""

    i = 0
    dictio = {}
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
                    dictio[str(sec_type)] = {str(seg) : (seg.x_xtra*1e-6, seg.y_xtra*1e-6, seg.z_xtra*1e-6)} #um to m
                else:
                    dictio[str(sec_type)][str(seg)] = (seg.x_xtra*1e-6, seg.y_xtra*1e-6, seg.z_xtra*1e-6) #um to m      

    return dictio


"""-----------------------------------------------------------------------------------WRITE REALNEURON (AND MODL)-----------------------------------------------------------------------------------"""
def read_mod(mech_folder, restrictions = None):
    """ read all .mod files in a directory and put them into a 2D-list
        restrictions: list that contains the mechanisms that need to be read
        :mech_folder: directory where the various NMODL mechanism files are stored
        :restrictions: a list of files that are the only ones that are considered in this function"""

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
    """read all .mod files in a directory and duplicate them in a seperate folder: eff
        :mech_folder: folder containing all the mechanisms .mod file
        :restrictions: selection of files that are duplicated"""

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
    """ read all all gbars from biophysics.hoc and put them into a dictionary
        :cell_folder: directory of the considered cell where the biophysics file is present
        :d_2_s: distance to the soma"""

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
    """ filters in the list of the model file texts the lines that contain a gating parameter
        :mod_files: list containing lists that contain all the lines of the NMODL .mod files that need to be translated
        :mod_names: list containing all the names of the different mechanisms"""

    l_alphas, l_betas, l_taus, l_infs = [], [], [], []
    hits = []
    for i,e in enumerate(mod_files):
        for j,f in enumerate(e):
            #print(test,f)
            l_alphas.append(([i,j],f,mod_names[i])) if re.search(tc.alpha_x_pattern,f) else None
            l_betas.append(([i,j],f,mod_names[i]))  if re.search(tc.beta_x_pattern,f)  else None
            l_taus.append(([i,j],f,mod_names[i]))   if re.search(tc.tau_x_pattern,f)   else None
            l_infs.append(([i,j],f,mod_names[i]))   if re.search(tc.x_inf_pattern,f)   else None
            hits.append(i) if  (re.search(tc.alpha_x_pattern,f) or re.search(tc.beta_x_pattern,f) or re.search(tc.tau_x_pattern,f) or re.search(tc.x_inf_pattern,f)) else None

    return l_alphas, l_betas, l_taus, l_infs, hits


def state_from_line(line,pattern):
    """ look if a gating state parameter is defined in the current/given line
        :line: line that can contain a gating state parameter
        :pattern: pattern that recognizes which kind of gating parameter state it is"""

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
    """ determine all the states from the 4 lists -> states names and descriptions
        :l_alphas: lines containing a formula with alpha(alpha_x_pattern)
        :l_betas: lines containing a formula with beta(beta_x_pattern)
        :l_taus: lines containing a formula with tau(tau_x_pattern)
        :l_infs: lines containing a formula with inf(x_inf_pattern)"""

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


def eq_hoc2pyt(equation):
    """ translate an equation or formula from hoc syntax to python syntax
        :equation: equation in hoc that needs to be translated to python"""

    equation = equation.replace('exp','np.exp').replace('^','**').replace('}','').replace('{','')
    equation = equation.strip() #remove spaces at beginning and end
    equation = equation.split(":")[0] #remove commentary in NMODL file
    equation = equation.lower() #put equations always in lowercase #POTENTIAL_RISK

    return equation


def gating_from_PROCEDURES_list(list_mod,mod_name): #,start_executing = 0):
    """ extract PROCEDURE block and save all equations in a list to retrieve alpha, beta, tau and inf from m and h
        :list_mod: list containing all the lines from the NMODL .mod file
        :mod_name: name of the mechanism"""

    #print(mod_name)
    celsius = tc.T_C #temperature in degrees Celsius
    proc_executing, param_executing = 0,0
    indents = 2

    equations_list = [4*indents*' '+'v = Vm', 4*indents*' '+f'celsius = {celsius}'] #uncomment this if v is a known parameter
    for i,e in enumerate(list_mod):
        estrip = e.strip()
    #for i,e in enumerate(list_mod):
        #stop when PROCEDURE BLOCK is finished
        if re.search('^}',e):
            proc_executing = 0
            param_executing = 0
        if proc_executing == 1:
            if re.search('if',e):    
                indents -= e.count("}")
                equations_list.append(4*indents*' '+eq_hoc2pyt(estrip)+":")
                indents += e.count("{")
            elif re.search('else',e):
                indents -= e.count("}")  
                equations_list.append(4*indents*' '+"else:")
                indents += e.count("{") 
            #elif e.search('if',e):
            elif re.search('=',e) and not(estrip.startswith(':')): #no if or else or fi or esle (end of)
                equation = eq_hoc2pyt(re.search(tc.equation_pattern,estrip).group())
                equations_list.append(4*indents*' '+equation)
            else:
                indents += e.count("{") - e.count("}") 
        elif param_executing == 1:
            if re.search('=',e) and not(estrip.startswith(':')):
                if re.findall("\(.*\)",e):
                    equations_list.append(4*indents*' '+eq_hoc2pyt(e.replace(re.findall("\(.*\)",e)[0],""))) #when evaluating the equations in PARAM block, the units are removed that are between brackets  
                else:
                    equations_list.append(4*indents*' '+eq_hoc2pyt(e)) #don't replace anything if units are not given
        #start when entering PROCEDURE block on the next iteration
        if re.search('PROCEDURE',e):
            proc_executing = 1
        #start when entering PARAMETER block on the next iteration
        if re.search('PARAMETER',e):
            param_executing = 1
    #maybe good idea to add a suffix to each variable -> NO
    #variables_dict[variable+'_'+mod_name.replace(".mod","")]
    for x in ['m','h']:
        if not (any([x+'alpha' in e for e in equations_list])) and (any([x+'tau' in e for e in equations_list])) and (any([x+'tau' in e for e in equations_list])):
            equations_list.append(4*indents*' '+f"{x+'alpha'} = {x+'inf'} / {x+'tau'} #only tau and inf provided in NMODL file")
        if not (any([x+'beta' in e for e in equations_list])) and (any([x+'tau' in e for e in equations_list])) and (any([x+'tau' in e for e in equations_list])):
            equations_list.append(4*indents*' '+f"{x+'beta'} = (1 - {x+'inf'}) / {x+'tau'} #only tau and inf provided in NMODL file")
    # if mod_name == 'Ca':
    #     print("final variables_dict: ",variables_dict)
    #     quit()
    
    return equations_list


def get_reversals(cell=None):
    """ extract the reversal voltages from the loaded Aberra cell in hoc
        this only works if a certain cell (sections) is loaded in
        :cell: dummy parameter that tells if a cell is loaded in or not"""
    if cell == None:
        return {'ek': [-85.0], 'ena': [50.0], 'eca': [132.4579341637009]}

    reversals = {"ek": [], "ehcn": [], "ena": [], "eca": []}

    for e in h.allsec():
        for rev, rev_vals in reversals.items():
            #print(e,rev)
            try:
                rev_val = eval(f"h.{e}.{rev}")
                if rev_val not in rev_vals:
                    if len(rev_vals) != 0:
                        print("multiple rev_vals for the same reversal potential?")
                    rev_vals.append(rev_val)
            except:
                pass

    return reversals


def currents_from_BREAKPOINT_list(list_mod,mod_name, gating_var):
    """ extract PROCEDURE block and save all equations to retrieve a list of the current from m and h
        :list_mod: a list containing all the lines of the NMODL .mod file
        :mod_name: the name of the mechanism
        :gating_var: gating parameters that are defined in the NMODL file"""

    celsius = tc.T_C #temperature in degrees Celsius
    break_executing, param_executing = 0,0
    indents = 2

    equations_list = [4*indents*' '+'v = Vm', 4*indents*' '+f'celsius = {celsius}'] #uncomment this if v is a known parameter
    reversals = get_reversals()
    for e,f in reversals.items():
        equations_list.append(4*indents*' '+f"{e} = {f[0]}")
    for e in gating_var:
        equations_list.append(4*indents*' '+f"{e} = {e}_{mod_name} #")

    for i,e in enumerate(list_mod):
        estrip = e.strip()
    #for i,e in enumerate(list_mod):
        #stop when PROCEDURE BLOCK is finished
        if re.search('^}',e):
            break_executing = 0
            param_executing = 0
        if break_executing == 1:
            if re.search('if',e):    
                indents -= e.count("}")
                equations_list.append(4*indents*' '+eq_hoc2pyt(estrip)+":")
                indents += e.count("{")
            elif re.search('else',e):
                indents -= e.count("}")  
                equations_list.append(4*indents*' '+"else:")
                indents += e.count("{") 
            #elif e.search('if',e):
            elif re.search('=',e) and not(estrip.startswith(':')): #no if or else or end of it: fi 
                equation = eq_hoc2pyt(re.search(tc.equation_pattern,estrip).group())
                equations_list.append(4*indents*' '+equation)
            else:
                indents += e.count("{") - e.count("}") 
        elif param_executing == 1:
            if re.search('=',e) and not(estrip.startswith(':')):
                if re.findall("\(.*\)",e):
                    equations_list.append(4*indents*' '+eq_hoc2pyt(e.replace(re.findall("\(.*\)",e)[0],""))) #when evaluating the equations in PARAM block, the units are removed that are between brackets  
                else:
                    equations_list.append(4*indents*' '+eq_hoc2pyt(e)) #don't replace anything if units are not given
        #start when entering PROCEDURE block on the next iteration
        if re.search('BREAKPOINT',e):
            break_executing = 1
        #start when entering PARAMETER block on the next iteration
        if re.search('PARAMETER',e):
            param_executing = 1
    equations_list.append(4*indents*' '+'vars='+str([eq.split('=')[0].strip()for eq in equations_list]))
    return equations_list


"""-----------------------------------------------------------------------------------WRITE MODL-----------------------------------------------------------------------------------"""
def one_to_multiline(root,file):
    """split a BLOCK defined a .modl file on one line over multiple lines {
    like this
    }
        :root:the directory where the file is stored
        :file:name of the file in the directory"""

    changed = 0
    with open(os.path.join(root,file), 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if re.search(tc.onelineBLOCK_pattern,line):
            part1, part2 = line.split("{")[0],line.split("{")[1].split("}")[0]
            lines[i] = f"{part1}{{\n\t{part2}\n}}"
            (line.split("{")[0],"\n",line.split("{")[1].split("}")[0])
            changed = 1
            print(f"{line} has been split over multiple lines in: {file}")
            #break
    if changed:
        with open(os.path.join(root,file), 'w') as newf: #.replace(".mod","_test.mod")
            newf.writelines(lines)

        return file #.replace(".mod","_test.mod")
    else:

        return file
    

def model_to_BLOCKS(flist):
    """ put the different BLOCK lines of a mechanism .modl file in a dictionary with each seperate block
        :flist: list containing all the lines of a file"""

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
    """ looks if a string is inside an element of a list(to search for certain keywords)
        :str: string is being looked for
        :list: list where string is searched"""

    return [e for e in list if str in e]


def str1_before_str2(str1,str2,flist):
    """ str1 NEEDS to be before str2 so interchange strings if str2 is before str1
        :str1: first string
        :str2: second string
        :flist: list containing all the lines of a file"""

    BLOCKdict = model_to_BLOCKS(flist)
    for e,f in BLOCKdict.items():
        list1 = str_in_element_in_list(str1,f) #local
        list2 = str_in_element_in_list(str2,f) #update
        if list1 and list2:
            loc1 = f.index(list1[0])
            loc2 = f.index(list2[0])
            if loc1 > loc2: #interchange them if str2 is before str1
                
                string1, string2 = f[loc1],f[loc2] #loc is location in BLOCK
                locat1 = flist.index(string1) #locat is location in whole flist
                locat2 = flist.index(string2) 
                flist[locat1], flist[locat2] = string2, string1
                return (flist, True)
    return (flist, False)

        # if any([str1 in g for g in f]) and any([str2 in g for g in f]):
        #     str1_index = f.index(str1)
        #     str2_index = f.index(str1)


def SUFFIX_Cm0(flist,Cm0):
    """ adds the capacitance value to the suffix to distinguish the different mechanisms in NEURON
        :flist: list containing the lines of the file that needs to be adapted
        :Cm0: membrane resting capacitance"""
    for i,e in enumerate(flist):
        if 'SUFFIX' in e:
            flist[i] = e.split('\n')[0]+Cm0+'\n' #put Cm0 value just before the newline
            break
    return flist


def eff_to_noteff(flist, Cm0):
    """ convert an effective version of a mech to a simple mechanism .mod file (this is the case for custom_pas)
        :flist: list containing the lines of the file that needs to be adapted
        :Cm0: membrane resting capacitance"""

    indices_to_remove = []
    for index, element in enumerate(flist): 
        if 'FUNCTION_TABLE' in element or 'INCLUDE' in element: #remove function tables and included files
            indices_to_remove.append(index)
    for index in sorted(indices_to_remove, reverse=True):
        del flist[index]

    for index, element in enumerate(flist): #replace the update line with a manual recast
        if 'update()' in element:
            flist[index] = f"Vm = v/{str(Cm0)} :Cm0 = 0.02 uF/cm2\n"
            flist.insert(index+1,'y = v\n')
    return flist
            

def replace_str(flist,to_repl,repl_with):
    """ replaces a certain string with another one in a given list
        :flist: list (possibly) containing the string to replace
        :to_repl: string that needs to be replaced with repl_with, this can also be a list
        :repl_with: string that needs to go in the place of to_repl"""
    if type(to_repl) == list:
        for i,e in enumerate(flist):
            for to,withh in zip(to_repl,repl_with):
                if to in e:
                    flist[i] = e.replace(to,withh)
                    break #if a replacement is done, the other strings in the list are ignored, only replace at most for 1 string in the list
        return flist
    for i,e in enumerate(flist):
        if to_repl in e:
            flist[i] = e.replace(to_repl,repl_with)
    return flist


"""-----------------------------------------------------------------------------------LUT MANIPULATION-----------------------------------------------------------------------------------"""
def merge_LUTdicts(dict1, dict2):
    """ merging 2 LUT dictionaries by combining their reference variable values and the corresponding tables
        :dict1: dictionary containing the first LUT
        :dict2: dicitonary containing the second LUT"""

    refs1, tables1 = dict1['refs'], dict1['tables']
    refs2, tables2 = dict2['refs'], dict2['tables']
    refs_merged, tables_merged = dict1['refs'].copy(), dict1['tables'].copy() #dict 1 and 2 contain the same refs (keys)

    for i,ref in enumerate(refs_merged): #iterate over the different reference variables/dimensions of the tables
        set1 = set(refs1[ref])
        set2 = set(refs2[ref])
        refs_merged[ref] = list(set1.union(set2))
        refs_merged[ref] = np.sort(refs_merged[ref])
        if len(refs1[ref]) != len(refs_merged[ref]) or len(refs2[ref]) != len(refs_merged[ref]): #refs1 and refs2 aren't fully overlapping so tables need to be combined
            for table in tables_merged:
                #tables_merged[table] = np.concatenate((tables1[table],tables2[table]),axis=i) #this assumes that refs1 and refs2 are disjoints sets
                                                                                               #and that the values of refs2 are all bigger than refs1
                for j, ref_val in enumerate(refs_merged[ref]):
                    index_val, whichrefs = (np.where(refs1[ref] == ref_val)[0][0], 1) if ref_val in refs1[ref] else (np.where(refs2[ref] == ref_val)[0][0], 2) if ref_val in refs2[ref] else None
                    #print(ref_val,index_val,whichrefs)
                    to_append = np.take(tables1[table],indices=range(index_val,index_val+1),axis=i) if whichrefs == 1 else np.take(tables2[table],indices=range(index_val,index_val+1),axis=i) if whichrefs == 2 else None
                    #print(to_append.shape)
                    new_table = to_append if j == 0 else np.concatenate((new_table,to_append),axis=i) 
                tables_merged[table] = new_table
                #check the new dimensions
                #print(tables_merged[table].shape, tables1[table].shape, tables2[table].shape)
        #print(len(set1),len(set2),len(refs_merged[ref]))

    merged_dict = {'refs': refs_merged, 'tables': tables_merged}

    return merged_dict


def merge_LUT(path, filename1, filename2):
    """ this reads in 2 LUT pickle files and writes out a LUT that combines the two
        :path: directory where both LUT files are stored
        :filename1: name of the first LUT file
        :filename2: name of the second LUT file"""

    pkldict1 = read_pickle(path,filename1)
    pkldict2 = read_pickle(path,filename2)
    merge_dict = merge_LUTdicts(pkldict1,pkldict2)
    print(filename1.split('.pkl'))
    filename1_split = filename1.split('.pkl')[0].split('_')
    filename2_split = filename2.split('.pkl')[0].split('_')
    merge_filename = ''
    for e,f in zip(filename1_split, filename2_split):
        if not(re.search('[a-zA-Z]',e) or re.search('[a-zA-Z]',f)): #if both filenames don't contain a letter, skip this part of the filename
                                                                    #relevant values need to include 
            continue
        if e == f:
            merge_filename += e
            merge_filename += '_'
            continue
        merge_filename += e + '_' + f + '_'   

    current_time = datetime.datetime.now()
    now = datetime.datetime.strftime(current_time,'%Y_%m_%d_%H_%M_%S')    
    merge_filename += now
    if merge_filename.endswith("_"):
        merge_filename = merge_filename[:-1]              
        
    load_pickle(merge_dict,path,merge_filename+'_merged.pkl')


def merge_LUTlist(path,filenamelist):
    """ this reads in a list of LUT pickle files and writes out a LUT that combines them all
        :path: directory where the various LUT files are stored
        :filenamelist: list of filenames of the LUTs that need to be merged"""

    pkldict_merged = read_pickle(path,filenamelist[0])
    merged_name = filenamelist[0]
    for e in filenamelist[1:]:

        pkldict2 = read_pickle(path,e)
        pkldict_merged = merge_LUTdicts(pkldict_merged,pkldict2)
        merged_split = merged_name.split('.pkl')[0].split('_')
        filename2_split = e.split('.pkl')[0].split('_')
        merge_filename = ''
        for e,f in zip(merged_split, filename2_split):
            if not(re.search('[a-zA-Z]',e) or re.search('[a-zA-Z]',f)): #if both filenames don't contain a letter, skip this part of the filename
                                                                        #relevant values need to include 
                continue
            if e == f:
                merge_filename += e
                merge_filename += '_'
                continue
            merge_filename += e + '_' + f + '_'   
        merged_name = merge_filename

    current_time = datetime.datetime.now()
    now = datetime.datetime.strftime(current_time,'%Y_%m_%d_%H_%M_%S')    
    merge_filename += now
    if merge_filename.endswith("_"):
        merge_filename = merge_filename[:-1]              
        
    load_pickle(pkldict_merged,path,merge_filename+'_merged.pkl')


def LUT_extend(filename):
    """ extends the LUT replacing the zero value by the edge values
        :filename: name of the LUT file that needs to be extended"""
    pkldict = read_pickle(filename)
    refs = pkldict['refs']
    for table in pkldict['tables'].values():
        shape = table.shape
        for i in range(shape[0]): #iterating over a
            for j in range(shape[1]): #iterating over f
                for k in range(shape[2]): #iterating over A

                    for m in range(shape[4]): #iterating over Cm0
                        for n in range(shape[5]): #iterating over fs
                            #print(table[i,j,k,:,m,n])
                            nonzero = np.nonzero(table[i,j,k,:,m,n])[0] #returns the indexes of nonzero values
                            first_value = table[i,j,k,nonzero[0],m,n]
                            last_value = table[i,j,k,nonzero[-1],m,n]
                            for l in range(shape[3]):
                                if l < nonzero[0]: #all values before the first nonzero value are put equal to that value
                                    table[i,j,k,l,m,n] = first_value
                                if l > nonzero[-1]: #same for all values above the last non-zero value
                                    table[i,j,k,l,m,n] = last_value
    load_pickle(pkldict,filename.replace('.pkl','_ext.pkl'))


def LUT_extend_1overtone(filename):
    """ extends the LUT replacing the zero value by the edge values
        :filename: name of the LUT file that needs to be extended"""
    pkldict = read_pickle(filename)
    refs = pkldict['refs']
    for table in pkldict['tables'].values():
        shape = table.shape
        for i in range(shape[0]): #iterating over a
            for j in range(shape[1]): #iterating over f
                for k in range(shape[2]): #iterating over A

                    for m in range(shape[4]): #iterating over AQ1
                        for n in range(shape[5]): #iterating over phiQ1
                            for o in range(shape[6]): #iterating over Cm0
                                for p in range(shape[7]): #iterating over fs                            
                            #print(table[i,j,k,:,m,n])
                                    nonzero = np.nonzero(table[i,j,k,:,m,n,o,p])[0] #returns the indexes of nonzero values
                                    first_value = table[i,j,k,nonzero[0],m,n,o,p]
                                    last_value = table[i,j,k,nonzero[-1],m,n,o,p]
                                    for l in range(shape[3]):
                                        if l < nonzero[0]: #all values before the first nonzero value are put equal to that value
                                            table[i,j,k,l,m,n,o,p] = first_value
                                        if l > nonzero[-1]: #same for all values above the last non-zero value
                                            table[i,j,k,l,m,n,o,p] = last_value
    load_pickle(pkldict,filename.replace('.pkl','_ext.pkl'))


def downsample_LUT(filename, down_factor=[1,1,2,2,1,1], load=True):
    """ downsamples an already calculated LUT
        :filename: location where the LUT that needs downsampling is located
        :down_factor: integer for each dimension showing the downsampling reduction factor"""
    
    pkldict = read_pickle(filename)
    refs = pkldict['refs']
    old_dims = [len(e) for e in refs.values()]
    #idea to only give specific dimension values
    # if new_dims:
    #     if len(new_dims) != 6:
    #         print('Wrong number of dimensions in downsampling')
    #         quit()
    if len(down_factor) != 6:
        raise ValueError('Wrong number of dimensions in downsampling')
    for i,ref in enumerate(pkldict['refs']): 
        pkldict['refs'][ref] = pkldict['refs'][ref][::down_factor[i]]
    for table in pkldict['tables']:
        pkldict['tables'][table] = pkldict['tables'][table][::down_factor[0], ::down_factor[1], ::down_factor[2],
                                        ::down_factor[3], ::down_factor[4], ::down_factor[5]]
    new_dims = [len(e) for e in refs.values()]
    print(f'dimension reduction: {old_dims} ---> {new_dims}')

    if load:
        load_pickle(pkldict,filename.replace('.pkl','_ds.pkl'))
    return pkldict


def upsample_LUT(filename, refs_up=None, method='linear', load=True):
    """ upsample an already calculated LUT
        :filename: location where the LUT that needs upsampling is located
        :up_factor: integer for each dimension showing the upsampling increase factor"""
    
    if not refs_up:
        a = np.array([16.0, 32.0, 64.0])*1e-9
        f = np.array([100., 500., 1e3, 2e3, 3e3])*1e3 #np.array([20., 100., 500., 1e3, 2e3, 3e3, 4e3])*1e3
        P_A = np.append(0,np.logspace(np.log10(0.1), np.log10(600), 50))*1e3
        Vm0min, Cm0max = -75, 0.02
        Qmin, Qmax = np.array([np.round(Vm0min - 35.0), 50.0]) * Cm0max * 1e-3
        DQ_LOOKUP = 1e-5
        Q = np.arange(Qmin, Qmax + DQ_LOOKUP, DQ_LOOKUP)
        Cm0 = np.array([0.01, 0.02])
        fs = np.array([0.75])
        refs_up = (a, f, P_A, Q, Cm0, fs)
    pkldict = read_pickle(filename)
    refs = pkldict['refs']
    old_dims = [len(e) for e in refs.values()]

    new_points = (a[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis, np.newaxis],
                    f[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis],
                    P_A[np.newaxis, np.newaxis, :, np.newaxis, np.newaxis, np.newaxis],
                    Q[np.newaxis, np.newaxis, np.newaxis, :, np.newaxis, np.newaxis],
                    Cm0[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :, np.newaxis],
                    fs[np.newaxis, np.newaxis, np.newaxis, np.newaxis, np.newaxis, :])
    # refs = tuple(e for e in refs.values())
    # print([e.shape for e in refs])
    # print('-'*150)
    # print([e.shape for e in refs_up])

    for table in pkldict['tables']:
        pkldict['tables'][table] = interp.interpn(tuple(e for e in refs.values()), pkldict['tables'][table],new_points, method=method)
        #print(f'done {table}')
    for i,ref in enumerate(pkldict['refs']): #this assumes that ref order is the same
        pkldict['refs'][ref] = refs_up[i]
    
    new_dims = [len(e) for e in refs_up]
    print(f'dimension increase: {old_dims} ---> {new_dims}')
    if load:
        load_pickle(pkldict,filename.replace('.pkl',f'_us_{method}.pkl'))
    return pkldict


def upsample_LUT2(filename, new_refs=None, method='linear', load=True):
    """ upsample an already calculated LUT
        :filename: location where the LUT that needs upsampling is located
        :up_factor: integer for each dimension showing the upsampling increase factor"""
    
    if not new_refs:
        a = np.array([16.0, 32.0, 64.0])*1e-9
        f = np.array([100., 500., 1e3, 2e3, 3e3])*1e3 #np.array([20., 100., 500., 1e3, 2e3, 3e3, 4e3])*1e3
        P_A = np.append(0,np.logspace(np.log10(0.1), np.log10(600), 50))*1e3
        Vm0min, Cm0max = -75, 0.02
        Qmin, Qmax = np.array([np.round(Vm0min - 35.0), 50.0]) * Cm0max * 1e-3
        DQ_LOOKUP = 1e-5
        Q = np.arange(Qmin, Qmax + DQ_LOOKUP, DQ_LOOKUP)
        Cm0 = np.array([0.01, 0.02])
        fs = np.array([0.75])
        new_refs = [a, f, P_A, Q, Cm0, fs]
    pkldict = read_pickle(filename)
    new_pkldict = copy.deepcopy(pkldict)
    old_refs = list((pkldict['refs']).values())
    old_dims = [len(e) for e in old_refs]
    new_dims = [len(e) for e in new_refs]
    print(f'old_dims: {old_dims}, new_dims: {new_dims}')


    new_points = (P_A[:, np.newaxis],
                    Q[np.newaxis, :])
    #ref_tuple = (refs['A'], refs['Q'], refs['Cm0'], refs['fs'])
    padding_vals = tuple((e-f,0) for e,f in zip(new_dims,old_dims))
    print(f"padding_vals: {padding_vals}")
    for table in pkldict['tables']:
        #print('before padding:', new_pkldict['tables'][table].shape)
        new_pkldict['tables'][table] = np.pad(pkldict['tables'][table],padding_vals)
        #print('after padding:', new_pkldict['tables'][table].shape)
        for i in range(len(pkldict['tables'][table])):
            for j in range(len(pkldict['tables'][table][i])):
                for m in range(len(pkldict['tables'][table][i,0,0,0])):
                    for n in range(len(pkldict['tables'][table][i,0,0,0,0])):
                        if method == 'akima':
                            akima_A = np.array([interp.Akima1DInterpolator(old_refs[3], row)(Q) for row in pkldict['tables'][table][i,j,:,:,m,n]])
                            akima = np.array([interp.Akima1DInterpolator(old_refs[2], col)(P_A) for col in akima_A.T]).T
                            new_pkldict['tables'][table][i,j,:,:,m,n] = akima
                        else:
                            new_pkldict['tables'][table][i,j,:,:,m,n] = interp.interpn(old_refs[2:4], pkldict['tables'][table][i,j,:,:,m,n], new_points, method=method)
        print(f'done {table}')
    for i,ref in enumerate(new_pkldict['refs']): #this assumes that ref order is the same
        new_pkldict['refs'][ref] = new_refs[i]
    
    print(f'dimension increase: {old_dims} ---> {new_dims}')
    if load:
        load_pickle(new_pkldict,filename.replace('.pkl',f'_us_{method}.pkl'))
    return new_pkldict


def random_LUT(filename):
    """ replaces the LUT with tables containing random values
        :filename: location where the LUT is stored"""
    pkldict = read_pickle(filename)
    for table in pkldict['tables']:
        pkldict['tables'][table] = np.random.rand(*pkldict['tables'][table].shape)
    load_pickle(pkldict,filename.replace('.pkl','_random.pkl'))


def lookup_LUT(filename, para = 'V', lookups = [32*1e-9, 500*1e3, 100*1e3, 50*1e-5, 0.01]):
    """ looks up the value in table given the reference values in each dimension
        :filename: location where the LUT is stored
        :para: lookup for this certain parameter table
        :lookups: the reference values for the LUT"""

    pkldict = read_pickle(filename)
    if len(lookups) != 5:
        raise LookupError("Wrong number of LUT values given")
    var_dict = {'a': lookups[0], 'f': lookups[1], 'A': lookups[2], 'Q': lookups[3], 'Cm0': lookups[4]}
    ind_list = [np.argmin(abs(pkldict['refs'][k]-v)) for (k,v) in var_dict.items()]
    vars = pkldict['refs']
    ref_list = [vars[e][f] for e,f in zip(['a','f','A','Q', 'Cm0'], ind_list)]
    print(f'lookup value for: {ref_list}')
    ind_list = f'[{ind_list[0]}, {ind_list[1]}, {ind_list[2]}, {ind_list[3]}, {ind_list[4]}, {0}]'
    value = eval(f"pkldict['tables']['{para}']{ind_list}")
    return value


def LUT_padding(filename):
    """"padding a LUT so it has the same dimensions in order to merge
        :filename: location where the unpadded LUT is located"""
    
    pkldict = read_pickle(filename)
    # for e in pkldict['refs']:
    #     print(e)
    pkldict['refs']['Q'] = np.pad(pkldict['refs']['Q'],(0,1))
    for e in pkldict['tables']:
        pkldict['tables'][e] = np.pad(pkldict['tables'][e],((0,0),(0,0),(0,0),(0,1),(0,0),(0,0)))
        #print(pkldict['tables'][e].shape)   
    load_pickle(pkldict,filename.replace('.pkl','_padded.pkl'))


def compare_LUT(filename1, filename2):
    """ compares 2 given LUT
        :filename1: location of the first LUT
        :filename2: location of the second LUT that needs to be compared with the first one"""
    pkldict1 = read_pickle(filename1)
    pkldict2 = read_pickle(filename2)
    V1 = pkldict1['tables']['V']#[0,0]
    V2 = pkldict2['tables']['V']
    if pkldict1['refs'].keys() != pkldict2['refs'].keys():
        # for e,f in zip(pkldict1['refs'].values(),pkldict2['refs'].values()):
        #     print(e,f,e-f)
        raise KeyError("reference values differ")
    refs = pkldict1['refs']
    #V1 = V1.reshape(V2.shape)
    print(f'shapes: {V1.shape}, {V2.shape}')
    print('[a, f, A, Q, Cm0, fs]')
    factors = [1e9, 1e-3, 1e-3, 1e5, 1e2,1e2]
    print('ABSOLUTE DIFFERENCE:')
    diff = abs(V1-V2)
    indexmin = np.unravel_index(np.argmin(diff),V2.shape)
    indexmax = np.unravel_index(np.argmax(diff),V2.shape)
    absmean, absmin, absmax = np.mean(diff), np.min(diff), np.max(diff)
    print(f'np.mean : {absmean}')
    print(f'np.min : {absmin} at : {[f"{e[f]*h:.2f}" for e,f,h in zip(refs.values(),indexmin,factors)]}, V1 : {V1[indexmin]}, V2 : {V2[indexmin]}')
    print(f'np.max : {absmax} at : {[f"{e[f]*h:.2f}" for e,f,h in zip(refs.values(),indexmax,factors)]}, V1 : {V1[indexmax]}, V2 : {V2[indexmax]}')
    print('-'*50)

    print('RELATIVE DIFFERENCE:')
    reldiff = abs(V1-V2)/abs(V2)
    indexrelmin = np.unravel_index(np.argmin(reldiff),V2.shape)
    indexrelmax = np.unravel_index(np.argmax(reldiff),V2.shape)
    relmean, relmin, relmax = np.mean(reldiff)*100, np.min(reldiff)*100, np.max(reldiff)*100
    relmeanp = f'{relmean:.2f}' if relmean < 100 else f'{relmean:.2e}'
    relminp = f'{relmin:.2f}' if relmin < 100 else f'{relmin:.2e}'
    relmaxp = f'{relmax:.2f}' if relmax < 100 else f'{relmax:.2e}'
    print(f'np.mean : {relmeanp} %') 
    print(f'np.min : {relminp} % at : {[f"{e[f]*h:.2f}" for e,f,h in zip(refs.values(),indexrelmin,factors)]}, V1 : {V1[indexrelmin]}, V2 : {V2[indexrelmin]}')
    print(f'np.max : {relmaxp} % at : {[f"{e[f]*h:.2f}" for e,f,h in zip(refs.values(),indexrelmax,factors)]}, V1 : {V1[indexrelmax]}, V2 : {V2[indexrelmax]}')
    print([e[f] for e,f in zip(refs.values(),indexmin)])
    print('\n')


def LUT_comparison(filename, factor = [1,1,2,2,1,1], method='linear', save=True):
    """ downsamples a LUT and upsamples it again to compare with the original LUT
        :filename: location where the (original) LUT is stored
        :method: interpolation method for the upsampling
        :save: the downsampled and usampled LUTs are saved in the same folder"""
    downsample_LUT(filename,factor)
    LUTdfile = filename.replace('.pkl','_ds.pkl')
    print('-'*100)
    upsample_LUT(LUTdfile, method=method)
    LUTufile = LUTdfile.replace('.pkl',f'_us_{method}.pkl')
    print('-'*100)
    compare_LUT(filename,LUTufile)
    if save:
        os.remove(LUTdfile)
        os.remove(LUTufile)


"""-----------------------------------------------------------------------------------PLOTTING-----------------------------------------------------------------------------------"""
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


def plt_LUTeff(variable,table,Q_refs,factor,unit,var_refs=None,plot=False,reduced_yrange=False, reduced_xrange=False):
    """ plots a certain gating parameter in function of the charge, possibility to provide multiple values for a tunable parameter
        :variable: the gating parameter that needs to be plotted
        :table: a table containing the y-values (for different values of a certain parameter)
        :Q_refs: x-values
        :factor: multiplication factor
        :unit: unit of the given parameter
        :var_refs: the parameter for which different plots is rendered given the different values
        :plot: showing the plot
        :reduced_yrange: reduction in y-range values
        :reduced_xrange: reduction in x-range values"""

    num_shades = len(var_refs) if type(var_refs) == np.ndarray else 1 #number of different colors for the different curves on 1 plot
    #actually it needs to be len(Arefs) but this results in a smaller range of colors which makes it more clear?
    step = 1 if num_shades < 10 else num_shades//10 #max 10 plots (sometimes it will be 11 due to rounding errors but OK)
    cmap = plt.cm.get_cmap('Wistia')
    # Generate colors across the colormap
    colors = [cmap(i / num_shades) for i in range(num_shades)]#.reverse()
    #colors.reverse()
    plt.clf()
    no_reduc = ['tcomp','V']
    for i,e in enumerate(table[::step]): #take only every fifth amplitude or other variable
        if num_shades > 1:
            plt.plot(Q_refs*1e5,e,label=f"{var_refs[::step][i]*factor:.1e} {unit}".replace('e+00','').replace('e+0','e').replace('e-0','e-'),color = colors[::step][i]) #plot Q with the (almost) 1D array of V (fs is still an extra dimension but this is no problem as fs only takes 1 value)
            plt.legend()
        else:
            plt.plot(Q_refs*1e5,table) #plot Q with the (almost) 1D array of V (fs is still an extra dimension but this is no problem as fs only takes 1 value)
                                        #plot the whole table as iterating over the elements will result in a reduced size of matrix
        plt.xlabel('$\mathrm{Q [nC/cm^2]}$')
        ylab = 'V [mV]' if variable == 'V' else 'C $\mathrm{[\dfrac{uF}{cm^2}]}$' if variable == 'C' else variable+ ' $\mathrm{[\dfrac{1}{ms}]}$' if ('alpha' in variable or 'beta' in variable) else variable
        plt.ylabel(f'{ylab}') #ylabel can be V, C or a gating parameter
        if reduced_yrange and (variable not in no_reduc): #only reduce the range for all gating variables and not for V or tcomp
            change_lower = plt.ylim()[0] < -10 #only change lower limit if values go below -10
            change_upper = plt.ylim()[1] > 10 #only change upper limit if values go above 10
            e_sorted = np.sort(e.reshape(-1)) #reshape(-1) just reshapes a multi-dimensional array to 1D
            e_limited = [e for e in e_sorted if e < 10 and e > -10] #filter values that go out of bounds
            if not e_limited: #if no value falls under the limited range, use the original array instead of no values
                e_limited = e
            # print(variable); print('before'); print(plt.ylim())
            plt.ylim(bottom = min(np.min(e_limited),-0.5)) if change_lower else None #lower value needs to be maximally -0.5
            plt.ylim(top = max(np.max(e_limited),0.5)) if change_upper else None #lower value needs to be minimally 0.5
        if reduced_xrange:
            plt.xlim(-100,50) #in order that the plots for different Cm0's have the same x-range (Q-range)
            # print('after'); print(plt.ylim()); print('\n\n')
        if plot:
            plt.show()  
        if num_shades < 2: #if only 1 plot needs to be shown, the 'whole' table is plotted and iterating over the different elements is useless as the elements lack a dimension
            break


def plt_LUTheat(variable,table,Q_refs, overtone_vars, factor,unit,var_refs=None,plot=False,reduced_yrange=False, reduced_xrange=False):
    """ plots a 2D 'heatmap' of a gating parameter in function of the charge Q and an overtone parameter
        :variable: the gating parameter that needs to be plotted
        :table: a table containing the y-values (for different values of a certain parameter)
        :Q_refs: x-values
        :factors: multiplication factor
        :units: unit of the given parameter
        :var_refs: the parameter for which different plots is rendered given the different values
        :plot: showing the plot
        :reduced_yrange: reduction in y-range values
        :reduced_xrange: reduction in x-range values"""

    # Plot the heatmap
    plt.clf()
    if table.ndim != 2:
        raise ValueError(f'Table has {table.ndim} dimensions but should have 2')
    plt.imshow(table, extent=(min(overtone_vars)*1e6, max(overtone_vars)*1e6, min(Q_refs)*1e5, max(Q_refs)*1e5,), origin='lower', cmap='viridis') 
    plt.colorbar(label=variable, orientation='horizontal')
    plt.xlabel('$\mathrm{AQ1 [10^{-1} nC/cm^2]}$')
    plt.ylabel('$\mathrm{Q [nC/cm^2]}$')
    if plot:
        plt.show()


def save_gatingplots(pkldict,foldername,reduced_yrange=True,reduced_xrange=False, a=32*1e-9, f=500*1e3, A=50*1e3, Cm0=0.01):
    """ plots the various gating parameters in function of the charge, it is possible to give multiple radius, frequencies, amplitudes or capacitances
        :pkldict: dictionary containing LUT
        :foldername: directory where the various plots will be stored
        :factor: multiplier to comply with given unit
        :unit: unit in which the parameter is expressed in
        :reduced_yrange: if a reduction in the gating parameter range is needed (most necessary when going to infinity)
        :reduced_xrange: to reduce the charge range when using different Cm0-values
        :a, f, A, Cm0: both specific values as the keyword 'all' can be given to specify which values need to be plugged in the independent parameters"""

    var_dict = {'a': a, 'f': f, 'A': A, 'Cm0': Cm0} #4 different variables can be sweeped over
    var_list = ['a', 'f', 'A', 'Cm0']
    all_list = np.where(np.array(list(var_dict.values()))=='all')[0] #looks which of the parameters has the 'all' keywoard
    if all_list: #only do this if there is a parameter that needs to be plotted for multiple values
        factor = tc.plotting_factors[var_list[all_list[0]]] #conversion factor from SI unit to a more conventional unit
        unit = tc.plotting_units[var_list[all_list[0]]] #the conventional unit
        title = var_list[all_list[0]]+': all' #the title starts with the parameter that contains the keywoard 'all' (for this parameter all values will be plotted)
    else:
        factor, unit = None, None
        title = ''
    for (k,v) in var_dict.items():
        if v != 'all':
            title += f", {k} = {v*tc.plotting_factors[k]:.1e} {tc.plotting_units[k]}" #each parameter is converted and then added with its respective conventional unit to the title
    title = title.replace('e+00','').replace('e+0','e').replace('e-0','e-') #shortening of scientific notation
    #all_list = [e=="all" for e in var_list.values()] #look which of the variables needs to be plotted over the whole range
    ind_list = [np.argmin(abs(pkldict['refs'][k]-v)) if v != 'all' else ':' for (k,v) in var_dict.items()] #where a single value is given, the array is indexed at the value that is closest to the given value
    print(f"used variables: {[pkldict['refs'][k][min] if min != ':' else 'all' for k, min in zip(var_dict,ind_list)]}")
    ind_list = ind_list[:3]+[":"]+[ind_list[3]]+[0] #a,f,A are indexed, Q needs to be plotted at the x-axis, Cm0 is indexed and fs can only have 1 value
    ind_list = [str(e) for e in ind_list] #convert the indexing list to strings for later on
    if len(all_list) > 1:
        raise ValueError('only 1 parameter can be plotted for all')
    Qrange = pkldict['refs']['Q'] #list containing the Q values on the x-axis
    var_refs = pkldict['refs'][var_list[all_list[0]]] if (len(all_list) > 0) else None #list containing all the values for the legend in the plot (each value results in a different curve)
    plt.figure(figsize=(8, 6)) #change figsize so all plots are shown properly and fit in the box

    #plot all calculated effective variables
    for key in pkldict['tables']: #iterate over the output gating parameters
        table = pkldict['tables'][key] #gating parameter matrix, dims:(a,f,A,Q,Cm,fs)
        table2 = eval(f"table[{(',').join(ind_list)}]") #index the table according to the given indexing list: 1 value for each parameter except all Q and a chosen variable
        if all_list:
            ind_colon = np.where(np.array(ind_list)==":")[0] #determine which variables keep all values (Q and a chosen variable)
            ind = ind_colon[ind_colon !=3][0] #except from Q (which has index 3), check the index of the variable that has all its values
            if ind > 3: #if the variable is further in the array than Q:
                table2 = np.swapaxes(table2,0,1) #switch the variable with the Q-axis so Q is plotted on the x-axis and the other one on the y-axis as desired
            table2 = np.reshape(table2,(-1,len(Qrange))) #reshape the table to dim(variable) x dim(Q) (this is partly done because of fs-dimension)
        plt_LUTeff(key,table2,Qrange,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #plotting
        plt.title(title)
        plt.savefig(f'figs/{foldername}/{key}.png')

    #plot the effective capacity
    table_V = eval(f"pkldict['tables']['V'][{(',').join(ind_list)}]") #index the table according to the calculated indexing list
    if all_list:
        if ind > 3: #if the variable is further in the array than Q:
            table_V = np.swapaxes(table_V,0,1) #switch the variable with the Q-axis so Q is plotted on the x-axis and the other one on the y-axis as desired
        table_V = np.reshape(table_V,(-1,len(Qrange))) #reshape the table to dim(variable) x dim(Q)
    #Q_table = np.tile(Qrange,np.prod(table_V.shape)//len(Qrange)).reshape(table_V.shape) #copy the Q array to a matrix with the same dimensions as V -> not needed, use broadcasting
    table_C = Qrange / table_V #Q = C x V
    plt_LUTeff('C',table_C,Qrange,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #plotting
    plt.title(title)
    plt.savefig(f'figs/{foldername}/C.png')  


def save_gatingplots_overtones(pkldict,foldername,reduced_yrange=True,reduced_xrange=False, a=32*1e-9, f=500*1e3, A=50*1e3, Cm0=0.01, overtones=1, overtone=1, heatmap=False):
    """ plots the various gating parameters in function of the charge, it is possible to give multiple radius, frequencies, amplitudes or capacitances
        :pkldict: dictionary containing LUT
        :foldername: directory where the various plots will be stored
        :factor: multiplier to comply with given unit
        :unit: unit in which the parameter is expressed in
        :reduced_yrange: if a reduction in the gating parameter range is needed (most necessary when going to infinity)
        :reduced_xrange: to reduce the charge range when using different Cm0-values
        :a, f, A, Cm0: both specific values as the keyword 'all' can be given to specify which values need to be plugged in the independent parameters
        :overtones: the number of input(=output in SONIC) overtones that the LUT contains
        :overtone: which overtone needs to be included in the plot"""

    var_dict = {'a': a, 'f': f, 'A': A, 'Cm0': Cm0} #4 different variables can be sweeped over
    var_list = ['a', 'f', 'A', 'Cm0']
    all_list = np.where(np.array(list(var_dict.values()))=='all')[0] #looks which of the parameters has the 'all' keywoard
    if all_list:
        factor = tc.plotting_factors[var_list[all_list[0]]] #conversion factor from SI unit to a more conventional unit
        unit = tc.plotting_units[var_list[all_list[0]]] #the conventional unit
        title = var_list[all_list[0]]+': all' #the title starts with the parameter that contains the keywoard 'all' (for this parameter all values will be plotted)
    else:
        factor, unit = None, None
        title = ''
    
    for (k,v) in var_dict.items():
        if v != 'all':
            title += f", {k} = {v*tc.plotting_factors[k]:.1e} {tc.plotting_units[k]}" #each parameter is converted and then added with its respective conventional unit to the title
    title = title.replace('e+00','').replace('e+0','e').replace('e-0','e-') #shortening of scientific notation
    #all_list = [e=="all" for e in var_list.values()] #look which of the variables needs to be plotted over the whole range
    ind_list = [np.argmin(abs(pkldict['refs'][k]-v)) if v != 'all' else ':' for (k,v) in var_dict.items()] #where a single value is given, the array is indexed at the value that is closest to the given value
    print(f"used variables: {[pkldict['refs'][k][min] if min != ':' else 'all' for k, min in zip(var_dict,ind_list)]}")
    ind_list = ind_list[:3]+[":"]+["0","0"]*(overtone-1)+[":","0"]+["0","0"]*(overtones-overtone)+[ind_list[3]]+[0] #a,f,A are indexed, Q needs to be plotted at the x-axis, Cm0 is indexed and fs can only have 1 value
                                                                                                                    #all overtones are indexed at 0,0 for both the amplitude and the phase, the given overtone has all amplitudes and the first phase
    ind_list = [str(e) for e in ind_list] #convert the indexing list to strings for later on
    if len(all_list) > 1:
        raise ValueError('only 1 parameter can be plotted for all')
    Qrange = pkldict['refs']['Q'] #list containing the Q values on the x-axis
    var_refs = pkldict['refs'][var_list[all_list[0]]] if (len(all_list) > 0) else None #list containing all the values for the legend in the plot (each value results in a different curve)
    overtone_amps = pkldict['refs'][f'AQ{overtone}'] #overtone amplitude array
    plt.figure(figsize=(8, 6)) #change figsize so all plots are shown properly and fit in the box

    #plot all calculated effective variables
    for key in pkldict['tables']: #iterate over the output gating parameters
        table = pkldict['tables'][key] #gating parameter matrix, dims:(a,f,A,Q,  A_Qx,phi_Qx  ,Cm,fs)
        table2 = eval(f"table[{(',').join(ind_list)}]") #index the table according to the given indexing list: 1 value for each parameter except all Q and a chosen variable
        if all_list:
            ind_colon = np.where(np.array(ind_list)==":")[0] #determine which variables keep all values (Q and a chosen variable)
            ind = ind_colon[(ind_colon<3) | (ind_colon>3+overtones)][-1] #except from Q and the overtones (which has index 3 up to 3 + 2*overtones), check the index of the variable that has all its values
            if ind > 3: #if the variable is further in the array than Q:
                table2 = table2.swapaxes(0,2).swapaxes(1,2) #reshuffle the axis so it has always the same order: variable, Q, AQx
            table2 = np.reshape(table2,(-1,len(Qrange),len(overtone_amps))) #reshape the table to dim(variable) x dim(Q) (this is partly done because of fs-dimension)
        if heatmap:
            plt_LUTheat(key,table2,Qrange,overtone_amps,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #2D heatmap
        else:
            plt_LUTeff(key,table2,Qrange,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #plotting
        plt.title(title)
        plt.savefig(f'figs/{foldername}/{key}.png')

    #plot the effective capacity
    table_V = eval(f"pkldict['tables']['V'][{(',').join(ind_list)}]") #index the table according to the calculated indexing list
    if all_list:
        if ind > 3: #if the variable is further in the array than Q:
            table_V = table_V.swapaxes(0,2).swapaxes(1,2) #reshuffle the axis so it has always the same order: variable, Q, AQx
        table_V = np.reshape(table_V,(-1,len(Qrange),len(overtone_amps))) #reshape the table to dim(variable) x dim(Q) x dim(AQx)
        table_C = (Qrange / table_V.swapaxes(2,1)).swapaxes(2,1) #axis are swapped and reswapped for proper broadcasting operation
    else:
        table_C = (Qrange / table_V.swapaxes(0,1)).swapaxes(0,1) #axis are swapped and reswapped for proper broadcasting operation
    if heatmap:
        plt_LUTheat('C',table_C,Qrange,overtone_amps,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #2D heatmap
    else:
        plt_LUTeff('C',table_C,Qrange,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #plotting
    plt.title(title)
    plt.savefig(f'figs/{foldername}/C.png')  


"""-----------------------------------------------------------------------------------DEPRECATED-----------------------------------------------------------------------------------"""
def formula_from_line(line, pattern, kinetic): #DEPRECATED?
    """split an equations into LHS and RHS"""
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


def formulas_from_lists(l_alphas, l_betas, l_taus, l_infs): #DEPRECATED?
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


def steadystates_from_gating_old(alphas, betas,taus,infs,dstates): #DEPRECATED
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


def steadystates_from_gating(states,gating_states_kinetics): #DEPRECATED?
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


def derstates_from_gating(states,gating_states_kinetics,x_dict): #DEPRECATED?
    """ derivative states functions are calculated based on the gating states kinetics and put into dictionary"""

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


def calc_eq(e,variables_dict,mod_name=None): #DEPRECATED
    """ parses an equation and puts the excecuted RHS into the variable of the LHS
        :e: equation that is parsed
        :variables_dict: all variables that are needed to execute the equation
        :mod_name: name of the mechanism"""

    LHS,RHS = e.split("=")
    variable = re.search(tc.var_pattern,LHS).group()
    formula = re.search(tc.math_pattern,RHS).group()
    equation = re.search(tc.equation_pattern,e).group()
    variable = eq_hoc2pyt(variable)
    formula = eq_hoc2pyt(formula)
    equation = eq_hoc2pyt(equation)
    #print(f'equation: {variable} : {formula} : {equation}')
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


def if_recursive(list_mod,if_statement,variables_dict,offset,mod_name=None): #DEPRECATED
    """ a recursive functions which handels executing if-statements encountered in a MODL file
        :list_mod: list containing the lines of a NMODL .mod file starting at an if case
        :if_statement: the value of the if_statement, indicating if the equation needs to be excecuted or not
        :variables_dict: dictionary containing all the executed variables (and also the ones needed to execute the equations)
        :offset: line offset indicating how many lines are being excecuted in if_recursive
        :mod_name: name of the mechanism"""

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


def gating_from_PROCEDURES_old(list_mod,mod_name,Vm): #,start_executing = 0): #DEPRECATED
    """ extract PROCEDURE block and execute all equations to retrieve alpha, beta, tau and inf from m and h
        :list_mod: list containing all the lines from a NMODL .mod file
        :mod_name: name of the mechanism
        :Vm: the voltage (v in neuron)"""

    #print(mod_name)
    celsius = tc.T_C #temperature in degrees Celsius
    proc_executing, param_executing = 0,0

    variables_dict = {'v': Vm, 'celsius' : celsius} #uncomment this if v is a known parameter
    # if mod_name == 'Ca':
    #     print(f"voltage = {variables_dict['v']}")
    i = 0
    while i < len (list_mod):
        e = list_mod[i]
        estrip = e.strip()
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
                if_statement = eval(eq_hoc2pyt(re.search(tc.if_pattern,e).group()),{**globals(),**variables_dict})
                #print(re.search(tc.if_pattern,e).group())
                # if mod_name == 'Ca':
                #     print('if_statement: ',if_statement)
                variables_dict,offset = if_recursive(list_mod[i+1:],if_statement,variables_dict,0,mod_name)
                # if mod_name == 'Ca':
                #     print(f"variables_dict : {variables_dict}, offset: {offset}")                    
                i += offset
                continue
            #elif?
            if re.search('=',e) and not(estrip.startswith(':')):
                #print(estrip)
                calc_eq(estrip,variables_dict,mod_name)
        elif param_executing == 1:
            if re.search('=',e) and not(estrip.startswith(':')):
                #print(estrip)
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


def currents_from_BREAKPOINT_old(list_mod,mod_name,Vm,x_dict,g_dict,location,start_executing = 0): #DEPRECATED
    """ extract PROCEDURE block and execute all equations to retrieve the current from m and h
        :list_mod: list containing the lines of the NMODL .mod file
        :mod_name: name of the mechanism
        :Vm: voltage (v in neuron)
        :x_dict: dictionary containing the gating parameters with their respective value of the mechanisms
        :g_dict:dictionary containing all the conductivities of the different mechanisms
        :location: keyword that reveals the location of the section which determines the exact conductivity
        :start_executing: parameter that determines if a line in the NMODL file needs to be executed or not"""

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


def plt_LUTeff_old(variable,table,Q_refs,A_refs,plot=False,reduced_yrange=False, reduced_xrange=False): #DEPRECATED
    """plot the given effective variable in function of the charge density for different pressure amplitudes
    this function assumes that the LUT is calculated for only 1 sonophore radius, 1 frequency and 1 sonophore coverage fraction"""

    num_shades = len(A_refs) #number of different colors for the different curves on 1 plot
    #actually it needs to be len(Arefs) but this results in a smaller range of colors which makes it more clear?
    step = 5

    cmap = plt.cm.get_cmap('Wistia')
    # Generate colors across the colormap
    colors = [cmap(i / num_shades) for i in range(num_shades)]#.reverse()
    #colors.reverse()
    #colors = plt.cm.Reds(np.linspace(0, 1, len(pkldict['refs']['Q'])))
    plt.clf()
    no_reduc = ['tcomp','V']
    for i,e in enumerate(table[::step]): #take only every fifth amplitude
        plt.plot(Q_refs*1e5,e,label=f"{A_refs[::step][i]*1e-3:.1e} kPa".replace('+0','').replace('-0',''),color = colors[::step][i]) #plot Q with the (almost) 1D array of V (fs is still an extra dimension but this is no problem as fs only takes 1 value)
        plt.legend()
        plt.xlabel('$\mathrm{Q [nC/cm^2]}$')
        ylab = 'V [mV]' if variable == 'V' else 'C $\mathrm{[\dfrac{uF}{cm^2}]}$' if variable == 'C' else variable+ ' $\mathrm{[\dfrac{1}{ms}]}$' if ('alpha' in variable or 'beta' in variable) else variable
        plt.ylabel(f'{ylab}')
        if reduced_yrange and (variable not in no_reduc): #only reduce the range for all gating variables and not for V or tcomp
            change_lower = plt.ylim()[0] < -10 #only change lower limit if values go below -10
            change_upper = plt.ylim()[1] > 10 #only change upper limit if values go above 10
            e_sorted = np.sort(e.reshape(-1)) #reshape(-1) just reshapes a multi-dimensional array to 1D
            e_limited = [e for e in e_sorted if e < 10 and e > -10] #filter values that go out of bounds
            if not e_limited: #if no value falls under the limited range, use the original array instead of no values
                e_limited = e
            # print(variable); print('before'); print(plt.ylim())
            plt.ylim(bottom = min(np.min(e_limited),-0.5)) if change_lower else None
            plt.ylim(top = max(np.max(e_limited),0.5)) if change_upper else None
        if reduced_xrange:
            plt.xlim(-100,50) #in order that the plots for different Cm0's have the same x-range (Q-range)
            # print('after'); print(plt.ylim()); print('\n\n')
        if plot:
            plt.show()
        

def save_gatingplots_old(pkldict,foldername,reduced_yrange=True,reduced_xrange=False,Cm0=None): #DEPRECATED
    """save all plots for the different LUTs, containing the effective variables of the gating kinetics
    (this function assumes that the LUT is calculated for only 1 sonophore radius, 1 frequency and 1 sonophore coverage fraction)"""

    Qrange,Arange = pkldict['refs']['Q'], pkldict['refs']['A']
    plt.figure(figsize=(8, 6)) #change figsize so all plots are shown properly and fit in the box
    #plot all calculated effective variables
    for key in pkldict['tables']:
        table = pkldict['tables'][key] #(a,f,A,Q,Cm,fs)
        table2 = table[0][0] #remove a and f
        if Cm0:
            ind = np.where(pkldict['refs']['Cm0']==Cm0)[0]
            if len(ind):
                ind = ind[0]
                table2 = np.moveaxis(table2,-2,0) #move Cm to the beginning
                table2 = table2[ind] #remove Cm
            else:
                print(f"Cm0 value: {Cm0} is not in LUT!\nPossible Cm0-values:{pkldict['refs']['Cm0']}")
                quit()
        plt_LUTeff_old(key,table2,Qrange,Arange,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange)
        Cm0_ext = str(Cm0) if Cm0 else '' #'_' + 
        plt.savefig(f'figs/{foldername}/{Cm0_ext}/{key}.png')

    #plot the effective capacity
    table_V = pkldict['tables']['V'][0][0]
    if Cm0:
        table_V = np.moveaxis(table_V,-2,0)
        table_V = table_V[ind]
    Q_table = np.tile(Qrange,np.prod(table_V.shape)//len(Qrange)).reshape(table_V.shape)
    table_C = Q_table / table_V
    plt_LUTeff_old('C',table_C,Qrange,Arange,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange)
    plt.savefig(f'figs/{foldername}/{Cm0_ext}/C.png')   