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
import csv
from scipy.signal import argrelextrema
import pprint
from scipy.optimize import minimize

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
    try:
        print(f'loaded pickle file with shape: {dictio["tables"]["V"].shape} in: {fpath}')
    except:
        pass
    

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


def read_csv(csv_file, print_time=0):

    with open(csv_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader: #iterate trough it as first row cannot be accessed that easily from the reader
            titles = row[1:] #first string in header is empty (row index)
            break
    sim_csv = np.genfromtxt(csv_file, delimiter=',',skip_header=1)
    sim_csv = sim_csv.swapaxes(1,0)[1:] #first column (now row after swapping) contains row index
    titles = [e.split(r"'")[1]+"_"+e.split(r"'")[3] if ('(' in e and ')' in e) else e for e in titles] #the gating parameters are stored as (x, mech), convert it to: x_mech
    factors =  [1e3 if e=='t' else 1 if e=='stimsate' else 1 if e=='Qm' else 1 if e=='Vm' else 1 if (e.startswith('i') or e.startswith('I')) else 1e2 if 'Cm' in e else 1 if ('_' in e) else 1 for e in titles]
    labels = ['t [ms]' if e=='t' else 'stimstate' if e=='stimstate' else 'Q [nC/cm2]' if e=='Qm' else 'V [mV]' if e=='Vm' else f'{e} [mA/m2]' if (e.startswith('i') or e.startswith('I')) else 'Cm [uF/cm2]' if 'Cm' in e else f'{e}' if ('_' in e) else '' for e in titles]

    if print_time:
        print(f"sim step: {sim_csv[0][2]*1e6} us, end sim: {sim_csv[0][-1]*1e3:.3f} ms")
    return sim_csv, titles, factors, labels


def create_plot_dict(titles, factors, labels, vars, variables, sim_csv):

    plot_dict = {}
    loc = 0
    for label, factor, array, title, var in zip(labels, factors, sim_csv, titles, vars):
        if 'nC/cm2' in label:
            Qrow = (array*factor)
        if 'mV' in label:
            Vrow = (array*factor)
        if var in variables:
            if var not in plot_dict.keys():
                plot_dict[var] = {'label': label, 'array': array, 'title': title, 'factor': factor}
                if var.startswith('m_') or var.startswith('h_'):
                        suffix = var.split('_')[-1] 
                        if 'm_'+ suffix in plot_dict.keys():
                            plot_dict['m_'+ suffix]['loc'] = loc
                        elif 'm_'+ suffix in variables:
                            plot_dict['m_'+ suffix] = {'loc': loc}
                        if 'h_'+ suffix in plot_dict.keys():
                            plot_dict['h_'+ suffix]['loc'] = loc
                        elif 'h_'+ suffix in variables:
                            plot_dict['h_'+ suffix] = {'loc': loc}
                        loc += 1
                elif (not var.startswith('i')) and (var != 't') and (not var.startswith('g')):
                    plot_dict[var]['loc'] = loc
                    loc += 1
            else:
                plot_dict[var].update({'label': label, 'array': array, 'title': title, 'factor': factor})
    if '' in plot_dict:
        del plot_dict['']
        loc -=1
    for e in plot_dict: #plot all the currents at the end
        if e.startswith('i'):
            plot_dict[e]['loc'] = loc

    if sum([e.startswith('g') for e in plot_dict]):
        loc +=1
        for e in plot_dict: #plot all the conductances at the end
            if e.startswith('g'):
                plot_dict[e]['loc'] = loc
    return plot_dict, Qrow, Vrow, loc


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


def add_dict_entry(dictio, value, keys):
    """adds an nested entry in a dictionary"""
    entry = dictio
    for i,key in enumerate(keys[:-1]):
        if key in entry.keys():
            entry = entry[key]
        else:
            for k in reversed(keys[i+1:]):
                value = {k:  value}
            entry[key] = value
            return dictio
    if keys[-1] in entry.keys():
        if entry[keys[-1]] != value:
            print(f'value overriden for {keys}: old value={entry[keys[-1]]}, new value={value}')
    entry[keys[-1]] = value
    return dictio


def txt_to_titration(txtfile,pklfile,save):
    """txt to a titration pkl file"""
    read = 1
    pkldict = read_pickle(pklfile)
    #pprint.pprint(pkldict)
    with open(txtfile) as tittxt:
        for i,line in enumerate(tittxt.readlines()):
                # if '---' in line:
                #     read = not read
                #     continue
                if read:
                    if len(line) < 2:
                        continue
                    line_splitted = line.split(':')
                    act = 1 if len(line_splitted) == 4 else 0 if len(line_splitted) == 2 else -1 #checks if the activation site is given
                    if act < 0:
                        print(f'wrong line: {line_splitted}')
                    if act:
                        values, remainder, section, time = line.split(':')
                        section, time = section.split(' ')[1], time.split(' ')[1]
                    else:
                        values, remainder = line.split(':')

                    if values[0] == '(' and values[-1] == ')':
                        values = values[1:-1]
                    else:
                        raise ValueError(f'Wrong values at {i}: {line}')
                    values = values.split(',')
                    #values: cell - freq - radius - coverage - DC - PRF
                    values = [float(f"{float(e):.2g}") for e in values] #round all values up to 2 decimals
                    values = [int(values[-1]), *values[:-1]] #cell can only be integer
                    if 'unexcitable' in remainder:
                        amp = -10*1e3
                        if act:
                            amp = (amp,section,float(time))
                        add_dict_entry(pkldict, amp, values)
                    else:
                        amp = float(remainder.split('kPa')[0])*1e3
                        if act:
                            amp = (amp,section,float(time))
                        add_dict_entry(pkldict, amp, values)

    #pprint.pprint(pkldict)
    if save:
        load_pickle(pkldict,pklfile)


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


def read_gbars(cell_folder, d_2_s):
    """ read all all gbars from biophysics.hoc and put them into a dictionary
        :cell_folder: directory of the considered cell where the biophysics file is present
        :d_2_s: distance to the soma"""

    #d_2_s = 0.001 #distance to soma is taken here as a constant but should be adapted to each section

    gbar_dict = {}
    #lines_list = []
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
        # if 'pas' in line and '=' in line:
        #     equation = eq_hoc2pyt(re.search(tc.equation_pattern,line.strip()).group())
        #     if 'pas' in gbar_dict.keys():
        #         gbar_dict['pas'] = np.append(gbar_dict['pas'],equation)
        #     else:
        #         gbar_dict['pas'] = [equation]
    return gbar_dict 


def read_biophysics(cell_folder, d_2_s):
    """ read all all gbars from biophysics.hoc and put them into a dictionary
        :cell_folder: directory of the considered cell where the biophysics file is present
        :d_2_s: distance to the soma"""

    #d_2_s = 0.001 #distance to soma is taken here as a constant but should be adapted to each section

    var_dict = {}
    #lines_list = []

    #read the biophysics.hoc file
    for root, _, files in os.walk(cell_folder):
        for file in files:
            if file.endswith("biophysics.hoc"):
                with open(os.path.join(root,file)) as f:
                    lines_list = f.readlines()

    #save the variables into a dictionary
    curr_loc = None
    for line in lines_list:
        search = re.search('forsec.*{',line) #search the start of a section type
        location = re.search("\$o1\.[a-zA-Z]*",line) #extract the section type
        if search and location:
            location = location.group().replace('$o1.','').lower()
            if curr_loc:
                raise TypeError(f'Already an assigned location: {curr_loc} so not possible to change it to: {location}')
            curr_loc = location #now we are inside {} brackets
            continue
        elif '}' in line:
            curr_loc = None #now we are oustide {} brackets
            continue
        elif curr_loc and '=' in line:
            equation = eq_hoc2pyt(line) #convert the equation from hoc to python
            var,expression = equation.split('=') #splot the equation in a LHS and RHS
            #if curr_loc == 'all': 
            #    var_dict[var] = eval(expression) #put the parameters directly into the dictionary if it is for all parameters (main reason is for the passive mechanism pas)
            #otherwise put it in a specific section type dictionary
            if curr_loc in var_dict: 
                var_dict[curr_loc][var.strip()] = eval(expression.strip())
            else:
                var_dict[curr_loc.strip()] = {var : eval(expression.strip())}
    return var_dict 


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
        if re.search('PROCEDURE',e): #to get the gating parameter formulas
            proc_executing = 1
        #start when entering PARAMETER block on the next iteration
        if re.search('PARAMETER',e): #to get the conductance
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
        #start when entering BREAKPOINT block on the next iteration
        if re.search('BREAKPOINT',e):
            break_executing = 1
        #start when entering PARAMETER block on the next iteration
        if re.search('PARAMETER',e):
            param_executing = 1
    variables = [eq.split('=')[0].strip()for eq in equations_list] #the variables of all the equation lines that will be written to the current method are separated
    #equations_list.append(4*indents*" "+"variables="+str(variables))
    currents = [e for e in variables if (e.startswith('i') or e.startswith('I'))] #from all the equations, there should be exactly 1 variable that represents the current
    if len(currents) != 1: #we expect 1 current to be calculated
        raise TypeError(f'Wrong amount of currents: {currents}')
    equations_list.append(4*indents*" "+"return "+currents[0]) #this current will be returned at the end of the method
    return equations_list


def conductances_from_BREAKPOINT_list(list_mod,mod_name, gating_var):
    """ extract PROCEDURE block and save all equations to retrieve a list of the conductance from m and h
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
        #start when entering BREAKPOINT block on the next iteration
        if re.search('BREAKPOINT',e):
            break_executing = 1
        #start when entering PARAMETER block on the next iteration
        if re.search('PARAMETER',e):
            param_executing = 1
    variables = [eq.split('=')[0].strip()for eq in equations_list] #the variables of all the equation lines that will be written to the current method are separated
    #equations_list.append(4*indents*" "+"variables="+str(variables))
    conductances = [e for e in variables if (e.startswith('g') and 'bar' not in e)] #from all the equations, there should be exactly 1 variable that represents the current
    if len(conductances) != 1: #we expect 1 current to be calculated
        print(f'No conductance found in {mod_name}.')
        equations_list.append(4*indents*" "+"return None")
    else:
        equations_list.append(4*indents*" "+"return "+conductances[0]) #this current will be returned at the end of the method
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
    if type(to_repl) == list: #to_repl and repl_with are both arrays containing a mapping
        for i,e in enumerate(flist): #iterate over the lines in the file
            for to,withh in zip(to_repl,repl_with): #iterate over the possible changes (mapping lists)
                if to in e: #if the line contains a 'to_change' string
                    e = e.replace(to,withh) #replace the string with the replacement
                    flist[i] = e.replace('g'+withh,'g'+to) #undo the ones that contain a conductance (conductance is the same for both mechanisms)
                    break #if a replacement is done, the other strings in the list are ignored, only replace at most for 1 string in the list
        return flist
    for i,e in enumerate(flist): #iterate over the lines in the file
        if to_repl in e: #if the string is in the line
            e = e.replace(to_repl,repl_with) #replace this line
            flist[i] = e.replace('g'+repl_with,'g'+to_repl) #undo conductances
    return flist


def FT_ov(flist,overtones):
    """ adds the extra arguments to the FUNCTION_TABLES
        :flist: list containing the lines of the file that needs to be adapted
        :overtones: the number of overtones"""
    extra_args = ''
    FUNCTION_TABLES = ''
    for ov in range(overtones):
        extra_args += f', Q{ov+1}(nC/cm2), phi{ov+1}(rad)'
        FUNCTION_TABLES += f'FUNCTION_TABLE A_V{ov+1}(A(kPa), Q(nC/cm2)) (mV)\nFUNCTION_TABLE phi_V{ov+1}(A(kPa), Q(nC/cm2)) (rad)\n'
    for i,e in enumerate(flist):
        if 'FUNCTION_TABLE' in e:
            flist[i] = e.replace('))',f'){extra_args})') 
            if flist[i+1].strip() == '':
                flist[i] += FUNCTION_TABLES.replace('))',f'){extra_args})') 
    return flist


def add_custom_pas(flist,overtones):
    """copied from write_modl_overtones.py because different for passive mechanism (custom_pas)"""
    overtone_ASSIGNED = ''
    overtone_NEURON = '    RANGE '
    for overtone in range(overtones):
        overtone_ASSIGNED += f'    q{overtone+1}  (nC/cm2)\n'
        overtone_ASSIGNED += f'    f{overtone+1}  (rad)\n'
        overtone_NEURON += f'q{overtone+1}, f{overtone+1}'

    block = None
    for i,e in enumerate(flist):
        if re.search(tc.block_init_pattern,e): #do determine if we are in a specific block or not
            block = re.search(tc.block_pattern,(re.search(tc.block_init_pattern,e).group(0))).group(0) #first look if we are in a BLOCK initiation line and then extract the actual block
        if block == "ASSIGNED" and flist[i+1].startswith('}'): #add extra lines at the end of the ASSIGNED block
            flist[i] = e + f'{overtone_ASSIGNED}'
        if block == "NEURON" and flist[i+1].startswith('}'): #add extra lines at the end of the NEURON block
            flist[i] = e + f'{overtone_NEURON}\n'
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


def lookup_LUT(filename, para = 'V', lookups = [32*1e-9, 500*1e3, 100*1e3, 40*1e-5, 0.75],Lemaire=False):
    """ looks up the value in table given the reference values in each dimension
        :filename: location where the LUT is stored
        :para: lookup for this certain parameter table
        :lookups: the reference values for the LUT
        :Lemaire: if a Lemaire-LUT is used"""

    pkldict = read_pickle(filename)
    #print(pkldict['refs'].keys())
    #print(pkldict['refs'].values())
    if len(lookups) != len(pkldict['refs'].keys()): #lookups bevat de waarde voor elke invoer parameter van de lookup
        raise LookupError(f"Wrong number of LUT values given: {len(lookups)} instead of {len(pkldict['refs'].keys())}")
    
    points = tuple(pkldict['refs'].values())
    interpolator = interp.RegularGridInterpolator(points,pkldict['tables'][para],method='linear')
    value = interpolator(lookups)

    #nearest neighbor own implementation
    # if Lemaire: 
    #     var_dict = {'a': lookups[0], 'f': lookups[1], 'A': lookups[2], 'Q': lookups[3], 'fs': lookups[4]}
    # else:
    #     var_dict = {'a': lookups[0], 'f': lookups[1], 'A': lookups[2], 'Q': lookups[3], 'Cm0': lookups[4]}
    # ind_list = [np.argmin(abs(pkldict['refs'][k]-v)) for (k,v) in var_dict.items()]
    # vars = pkldict['refs']
    # if Lemaire:
    #     ref_list = [vars[e][f] for e,f in zip(['a','f','A','Q', 'fs'], ind_list)]
    # else:
    #     ref_list = [vars[e][f] for e,f in zip(['a','f','A','Q', 'Cm0'], ind_list)]
    # print(f'lookup value for: {ref_list}')
    # if Lemaire:
    #     ind_list = f'[{ind_list[0]}, {ind_list[1]}, {ind_list[2]}, {ind_list[3]}, {ind_list[4]}]'
    # else:
    #     ind_list = f'[{ind_list[0]}, {ind_list[1]}, {ind_list[2]}, {ind_list[3]}, {ind_list[4]}, {0}]'
    # value = eval(f"pkldict['tables']['{para}']{ind_list}")
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


def compare_LUT(filename1, filename2, nvar=5, violin=False):
    """ compares 2 given LUT
        :filename1: location of the first LUT -> original (actual values)
        :filename2: location of the second LUT that needs to be compared with the first one
        :nvar: number of variables that are considered"""
    
    ABSOLUTE, RELATIVE = 0, 1
    ivar = 0
    pkldict1 = read_pickle(filename1)
    pkldict2 = read_pickle(filename2)
    if pkldict1['refs'].keys() != pkldict2['refs'].keys():
        # for e,f in zip(pkldict1['refs'].values(),pkldict2['refs'].values()):
        #     print(e,f,e-f)
        raise KeyError("reference values differ")
    print(f"shapes: {pkldict1['tables']['V'].shape}, {pkldict2['tables']['V'].shape}")
    print('     [  a,      f,       A,      Q,      Cm0,     fs]')
    refs = pkldict1['refs']
    factors = [1e9, 1e-3, 1e-3, 1e5, 1e2,1e2]
    print(f"min: {[f'{e[0]*f:.1f}' for e,f in zip(refs.values(), factors)]}")
    print(f"max: {[f'{e[-1]*f:.1f}' for e,f in zip(refs.values(), factors)]}")
    for var, V1, V2 in zip(pkldict1['tables'], pkldict1['tables'].values(), pkldict2['tables'].values()):
        #V1 = V1.reshape(V2.shape)
        print(var)
        if ABSOLUTE:
            print('ABSOLUTE DIFFERENCE:')
            diff = abs(V1-V2)
            indexmin = np.unravel_index(np.argmin(diff),V2.shape)
            indexmax = np.unravel_index(np.argmax(diff),V2.shape)
            absmean, absmin, absmax = np.mean(diff), np.min(diff), np.max(diff)
            #print(f'np.mean : {absmean}')
            #print(f'np.min : {absmin} at : {[f"{e[f]*h:.2f}" for e,f,h in zip(refs.values(),indexmin,factors)]}, V1 : {V1[indexmin]}, V2 : {V2[indexmin]}')
            print(f'np.max : {absmax} at : {[f"{e[f]*h:.2f}" for e,f,h in zip(refs.values(),indexmax,factors)]}, V1 : {V1[indexmax]}, V2 : {V2[indexmax]}')
            print('-'*50)

        if RELATIVE:
            print('RELATIVE DIFFERENCE:')
            reldiff = abs(V1-V2)/abs(V1)
            indexrelmin = np.unravel_index(np.argmin(reldiff),V2.shape)
            indexrelmax = np.unravel_index(np.argmax(reldiff),V2.shape)
            relmean, relmed, relmin, relmax = np.mean(reldiff)*100, np.median(reldiff)*100, np.min(reldiff)*100, np.max(reldiff)*100
            percentile = 75
            relpercentile = np.percentile(reldiff,percentile)*100
            relmeanp = f'{relmean:.2f}' if relmean < 1e3 else f'{relmean:.2e}'
            relmedp = f'{relmed:.2f}' if relmed < 1e3 else f'{relmed:.2e}'
            relperp = f'{relpercentile:.2f}' if relpercentile < 1e3 else f'{relpercentile:.2e}'
            relminp = f'{relmin:.2f}' if relmin < 1e3 else f'{relmin:.2e}'
            relmaxp = f'{relmax:.2f}' if relmax < 1e3 else f'{relmax:.2e}'
            #print(f'np.mean : {relmeanp} %') 
            print(f'np.median : {relmedp} %')
            print(f'np.percentile({percentile}) : {relperp} %')
            #print(f'np.min : {relminp} % at : {[f"{e[f]*h:.2f}" for e,f,h in zip(refs.values(),indexrelmin,factors)]}, V1 : {V1[indexrelmin]}, V2 : {V2[indexrelmin]}')
            print(f'np.max : {relmaxp} % at : {[f"{e} = {refs[e][f]*tc.all_factors[e]:.2f}" for e,f in zip(refs.keys(),indexrelmax)]}, {var}1 : {V1[indexrelmax]}, {var}2 : {V2[indexrelmax]}')
            if violin:
                plt.violinplot(reldiff.reshape(-1)*100)
                plt.show()
        print('\n')
        ivar += 1
        if ivar == nvar:
            break


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


def LUT_to_LUT2(filename, remove_zeros=False):
    """ converts a LUT gating parameters to a LUT where the values for Cm0=0.01 and 0.02 are separated into different parameters
        :filename: location of the original LUT
        :remove_zeros: remove the uncalculated 0-values in the 0.01-variant"""
    
    pkldict = read_pickle(filename)
    tables = list(pkldict['tables'].keys())
    for table in tables:
        pkldict['tables'][table+'2'] = pkldict['tables'][table][:,:,:,:,1,:]
        pkldict['tables'][table] = pkldict['tables'][table][:,:,:,:,0,:]
        if remove_zeros:
            #print(pkldict['tables'][table].shape)
            if table == 'tcomp': #this table doesn't contain zeros but needs to be reduced to comply with the shape
                pkldict['tables'][table] = pkldict['tables'][table][:,:,:,nonzeros]
                continue
            nonzeros = np.nonzero(pkldict['tables'][table][0,0,0,:,0])[0] # array is given inside a tuple hence the slicing
            pkldict['tables'][table] = pkldict['tables'][table][:,:,:,nonzeros] #remove zeros in the Q-axis
            nonzeros2 = np.nonzero(pkldict['tables'][table+'2'][0,0,0,:,0])[0] # array is given inside a tuple hence the slicing
            pkldict['tables'][table+'2'] = pkldict['tables'][table+'2'][:,:,:,nonzeros2] #remove zeros in the Q-axis
            #print(pkldict['tables'][table].shape)
    del pkldict['refs']['Cm0'] #delete reference key as this is not a dimension anymore
    if remove_zeros:
        pkldict['tables']['Q_ext'] = pkldict['refs']['Q'][nonzeros2]
        pkldict['refs']['Q'] = pkldict['refs']['Q'][nonzeros]
    #print(tables)
    load_pickle(pkldict,filename.replace('.pkl','_LUT2.pkl'))


def LUT_to_LUT2_1overtone(filename, remove_zeros=False):
    """ converts a LUT gating parameters to a LUT where the values for Cm0=0.01 and 0.02 are separated into different parameters
        :filename: location of the original LUT
        :remove_zeros: remove the uncalculated 0-values in the 0.01-variant"""
    
    i_vars = [0, 0, 11, 0, 0, 0] #use amplitude different than zero, otherwise zero values can be present in the LUT that are not dummy 0's
    a, f, A, A1, phi1, fs = i_vars
    pkldict = read_pickle(filename)
    tables = list(pkldict['tables'].keys())
    # for e,f in zip(pkldict['refs'],i_vars):
    #     print(e, pkldict['refs'][e][f])
    for table in tables:
        pkldict['tables'][table+'2'] = pkldict['tables'][table][:,:,:,:,:,:,1,:]
        pkldict['tables'][table] = pkldict['tables'][table][:,:,:,:,:,:,0,:]
        if remove_zeros:
            #print(pkldict['tables'][table].shape)
            if table == 'tcomp': #this table doesn't contain zeros but needs to be reduced to comply with the shape
                pkldict['tables'][table] = pkldict['tables'][table][:,:,:,nonzeros,:,:]
                continue
            nonzeros = np.nonzero(pkldict['tables'][table][a,f,A,:,A1,phi1,fs])[0] 
            #print(nonzeros)
            pkldict['tables'][table] = pkldict['tables'][table][:,:,:,nonzeros,:,:,:] #remove zeros in the Q-axis
            nonzeros2 = np.nonzero(pkldict['tables'][table+'2'][a,f,A,:,A1,phi1,fs])[0]
            pkldict['tables'][table+'2'] = pkldict['tables'][table+'2'][:,:,:,nonzeros2,:,:,:] #remove zeros in the Q-axis
            #print(pkldict['tables'][table+'2'].shape)
    del pkldict['refs']['Cm0'] #delete reference key as this is not a dimension anymore
    if remove_zeros:
        pkldict['tables']['Q_ext'] = pkldict['refs']['Q'][nonzeros2]
        pkldict['refs']['Q'] = pkldict['refs']['Q'][nonzeros]
    #print(tables)
    load_pickle(pkldict,filename.replace('.pkl','_LUT2.pkl'))


def LUTov_to_LUT(filename):
    """ remove the overtone variables of a LUT to compare with the original LUT with no overtones
        :filename: location of the original LUT"""
    
    pkldict = read_pickle(filename)
    refs = pkldict['refs']
    tables = pkldict['tables']
    ov = 1
    while len(refs) != 5:
        del refs['AQ'+str(ov)]
        del refs['phiQ'+str(ov)]
        for table in tables:
            if len(tables[table].shape) > 5:
                tables[table] = tables[table][:,:,:,:,0,0]
        ov += 1
    load_pickle(pkldict,filename.replace('.pkl','_0ov.pkl'))


def minimize_ov(filename_0ov,filename_1ov):
    """looks up for which overtone values the lookup value is the closest to the non-overtone lookup table: this should be for the null vector
        :filename_0ov: the non-overtone lookup
        :filename_1ov: the lookup table with overtones"""
    
    # lkp0 = read_pickle(filename_0ov)
    # lkp1 = read_pickle(filename_1ov)
    # for e,f in zip(lkp0['refs'].values(),lkp1['refs'].values()):
    #     print(e,f)
    def lookup(inp_vars):
        q1,f1 = inp_vars[0], inp_vars[1]
        val0 = lookup_LUT(filename_0ov,lookups=[32*1e-9, 500*1e3, 600*1e3, -102*1e-5, 0.75])
        val1 = lookup_LUT(filename_1ov,lookups=[32*1e-9, 500*1e3, 600*1e3, -102*1e-5, q1, f1, 0.75])
        #print(val,val2)
        values =  abs((val1-val0)[0])
        return values
    sol = minimize(lookup, [0, 0],bounds=([0,0.001],[0,5.026]))
    print(sol.x)
    val0 = lookup_LUT(filename_0ov,lookups=[32*1e-9, 500*1e3, 600*1e3, -102*1e-5, 0.75])
    val1 = lookup_LUT(filename_1ov,lookups=[32*1e-9, 500*1e3, 600*1e3, -102*1e-5, 0,0, 0.75])
    print(val0,val1)

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


def plt_LUTeff(variable,table,Q_refs,factor,unit,var_refs=None,plot=False,reduced_yrange=False, reduced_xrange=False, label = False):
    """ plots a certain gating parameter in function of the charge, possibility to provide multiple values for a tunable parameter
        :variable: the gating parameter that needs to be plotted
        :table: a table containing the y-values (for different values of a certain parameter)
        :Q_refs: x-values
        :factor: multiplication factor
        :unit: unit of the given parameter
        :var_refs: the parameter for which different plots is rendered given the different values
        :plot: showing the plot
        :reduced_yrange: reduction in y-range values
        :reduced_xrange: reduction in x-range values
        :label: specific label can be added if there is no var_refs given"""

    num_shades = len(var_refs) if type(var_refs) == np.ndarray else 1 #number of different colors for the different curves on 1 plot
    #actually it needs to be len(Arefs) but this results in a smaller range of colors which makes it more clear?
    step = 1 if num_shades < 10 else num_shades//10 #max 10 plots (sometimes it will be 11 due to rounding errors but OK)
    cmap = plt.cm.get_cmap('Wistia')
    # Generate colors across the colormap
    colors = [cmap(i / num_shades) for i in range(num_shades)]#.reverse()
    #colors.reverse()
    #plt.clf()
    no_reduc = ['tcomp','V','V2']
    for i,e in enumerate(table[::step]): #take only every fifth amplitude or other variable
        mul = 1 if (variable == 'V' or variable == 'V2') else 1e5 if (variable == 'C' or variable == 'C2') else 1 if ('alpha' in variable or 'beta' in variable) else 1 #factor is 1 for alpha and beta if in 1/ms, 1e3 if 1/s 
        if num_shades > 1:
            plt.plot(Q_refs*1e5,e*mul,label=f"{var_refs[::step][i]*factor:.1e} {unit}".replace('e+00','').replace('e+0','e').replace('e-0','e-'),color = colors[::step][i]) #plot Q with the (almost) 1D array of V (fs is still an extra dimension but this is no problem as fs only takes 1 value)
            plt.legend()
        elif label:
            if 'beta' in variable:
                if 'betah' in variable:
                    plt.plot(Q_refs*1e5,table*mul,'r--',label=variable)
                else:
                    plt.plot(Q_refs*1e5,table*mul,'y--',label=variable)
            else:
                if 'alpham' in variable:
                    plt.plot(Q_refs*1e5,table*mul,'y-',label=variable)
                else:
                    plt.plot(Q_refs*1e5,table*mul,'r-',label=variable)
        else:
            plt.plot(Q_refs*1e5,table*mul) #plot Q with the (almost) 1D array of V (fs is still an extra dimension but this is no problem as fs only takes 1 value)
                                        #plot the whole table as iterating over the elements will result in a reduced size of matrix
        plt.xlabel('$\mathrm{Q [nC/cm^2]}$')
        if 'alpha' in variable:
            var, suf = r'\alpha', variable.split('_')[0][-1]+", "+variable.split('_')[1] 
        elif 'beta' in variable:
            var, suf = r'\alpha', variable.split('_')[0][-1]+", "+variable.split('_')[1]
        ylab = '${V_m}^*$ [mV]' if variable == 'V' else '${C_m}^*$ $\mathrm{[\dfrac{uF}{cm^2}]}$' if variable == 'C' else f'${{{var}}}_{{{suf}}}^{{*}}$' + ' $\mathrm{[\dfrac{1}{ms}]}$' if ('alpha' in variable or 'beta' in variable) else variable
        plt.ylabel(f'{ylab}') #ylabel can be V, C or a gating parameter
        if reduced_yrange and (variable not in no_reduc): #only reduce the range for all gating variables and not for V or tcomp
            lower_lim, upper_lim = -20, 20
            change_lower = plt.ylim()[0] < lower_lim #only change lower limit if values go below -10
            change_upper = plt.ylim()[1] > upper_lim #only change upper limit if values go above 10
            e_sorted = np.sort(e.reshape(-1)) #reshape(-1) just reshapes a multi-dimensional array to 1D
            e_limited = [e for e in e_sorted if e < upper_lim and e > lower_lim] #filter values that go out of bounds
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

    var_dict = {'a': a, 'f': f, 'A': A, 'Cm0': Cm0} if Cm0 else {'a': a, 'f': f, 'A': A} #4 different variables can be sweeped over
    var_list = ['a', 'f', 'A', 'Cm0'] if Cm0 else ['a', 'f', 'A']
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
    ind_list = ind_list[:3]+[":"] + list(ind_list[3:]) + [0] if Cm0 else ind_list[:3]+[":"] + [0] #a,f,A are indexed, Q needs to be plotted at the x-axis, Cm0 is indexed, and fs can only have 1 value
    #only add Cm0 if LUT contains Cm0 dimension
    ind_list = [str(e) for e in ind_list] #convert the indexing list to strings for later on
    if len(all_list) > 1:
        raise ValueError('only 1 parameter can be plotted for all')
    Qrange = pkldict['refs']['Q'] #list containing the Q values on the x-axis
    if 'Q_ext' in pkldict['tables']:
        Q_ext = pkldict['tables']['Q_ext'] #define Q_ext if it is defined in the LUT
        del pkldict['tables']['Q_ext'] #remove it from the LUT as it is not needed to plot
    var_refs = pkldict['refs'][var_list[all_list[0]]] if (len(all_list) > 0) else None #list containing all the values for the legend in the plot (each value results in a different curve)
    plt.figure(figsize=(8, 6)) #change figsize so all plots are shown properly and fit in the box

    #plot all calculated effective variables
    for key in pkldict['tables']: #iterate over the output gating parameters
        plt.clf()
        if 'tcomp' in key:
            continue #don't plot the computation time tcomp as these plots are not important
        table = pkldict['tables'][key] #gating parameter matrix, dims:(a,f,A,Q,Cm,fs)
        Q_array = Qrange if len(Qrange) in table.shape else Q_ext #if the length of Qrange is the same as on of the dimensions in the table, use this array, otherwise use the extended Q array
        table2 = eval(f"table[{(',').join(ind_list)}]") #index the table according to the given indexing list: 1 value for each parameter except all Q and a chosen variable
        if all_list:
            ind_colon = np.where(np.array(ind_list)==":")[0] #determine which variables keep all values (Q and a chosen variable)
            ind = ind_colon[ind_colon !=3][0] #except from Q (which has index 3), check the index of the variable that has all its values
            if ind > 3: #if the variable is further in the array than Q:
                table2 = np.swapaxes(table2,0,1) #switch the variable with the Q-axis so Q is plotted on the x-axis and the other one on the y-axis as desired
            table2 = np.reshape(table2,(-1,len(Q_array))) #reshape the table to dim(variable) x dim(Q) (this is partly done because of fs-dimension)
        plt_LUTeff(key,table2,Q_array,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #plotting
        #plt.title(title)
        plt.savefig(f'figs/{foldername}/{key}.png')

    #plot the effective capacity
    table_V = eval(f"pkldict['tables']['V'][{(',').join(ind_list)}]") #index the table according to the calculated indexing list
    #print(pkldict['tables'].keys())
    t_V2 = 0
    if t_V2:
        table_V2 = eval(f"pkldict['tables']['V2'][{(',').join(ind_list)}]")
    if all_list:
        if ind > 3: #if the variable is further in the array than Q:
            table_V = np.swapaxes(table_V,0,1) #switch the variable with the Q-axis so Q is plotted on the x-axis and the other one on the y-axis as desired
        table_V = np.reshape(table_V,(-1,len(Qrange))) #reshape the table to dim(variable) x dim(Q)
    #Q_table = np.tile(Qrange,np.prod(table_V.shape)//len(Qrange)).reshape(table_V.shape) #copy the Q array to a matrix with the same dimensions as V -> not needed, use broadcasting
    table_C = Qrange / table_V#[:,0] #Q = C x V <=> C = Q / V: C/ m2 / mV = 1 / 1e-3 * C/m2/V = 1e3 F/m2 = 1e5 uF/m2
    if t_V2:
        table_C2 = Q_ext / table_V2
    plt.clf()
    plt_LUTeff('C',table_C,Qrange,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #plotting
    plt.title(title)
    #plt.show()
    plt.savefig(f'figs/{foldername}/C.png')  
    plt.clf()
    if t_V2:
        plt_LUTeff('C2',table_C2,Q_ext,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange) #plotting
        plt.title(title)
        plt.savefig(f'figs/{foldername}/C2.png') 


def save_gatingplots_group(pkldict,foldername,reduced_yrange=True,reduced_xrange=False, a=32*1e-9, f=500*1e3, A=50*1e3, Cm0=0.01):
    """ plots the various gating parameters, where they are grouped based on the mechanism, in function of the charge, it is possible to give multiple radius, frequencies, amplitudes or capacitances
        :pkldict: dictionary containing LUT
        :foldername: directory where the various plots will be stored
        :factor: multiplier to comply with given unit
        :unit: unit in which the parameter is expressed in
        :reduced_yrange: if a reduction in the gating parameter range is needed (most necessary when going to infinity)
        :reduced_xrange: to reduce the charge range when using different Cm0-values
        :a, f, A, Cm0: both specific values as the keyword 'all' can be given to specify which values need to be plugged in the independent parameters"""

    var_dict = {'a': a, 'f': f, 'A': A, 'Cm0': Cm0} if Cm0 else {'a': a, 'f': f, 'A': A} #4 different variables can be sweeped over
    var_list = ['a', 'f', 'A', 'Cm0'] if Cm0 else ['a', 'f', 'A']
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
    ind_list = ind_list[:3]+[":"] + list(ind_list[3:]) + [0] if Cm0 else ind_list[:3]+[":"] + [0] #a,f,A are indexed, Q needs to be plotted at the x-axis, Cm0 is indexed, and fs can only have 1 value
    #only add Cm0 if LUT contains Cm0 dimension
    ind_list = [str(e) for e in ind_list] #convert the indexing list to strings for later on
    if len(all_list) > 1:
        raise ValueError('only 1 parameter can be plotted for all')
    Qrange = pkldict['refs']['Q'] #list containing the Q values on the x-axis
    if 'Q_ext' in pkldict['tables']:
        Q_ext = pkldict['tables']['Q_ext'] #define Q_ext if it is defined in the LUT
        del pkldict['tables']['Q_ext'] #remove it from the LUT as it is not needed to plot
    var_refs = pkldict['refs'][var_list[all_list[0]]] if (len(all_list) > 0) else None #list containing all the values for the legend in the plot (each value results in a different curve)
    plt.figure(figsize=(8, 6)) #change figsize so all plots are shown properly and fit in the box

    group_dict = {}
    for key in pkldict['tables']:
        if '_' in key:
            pre,suf = key.split('_')
            if suf in group_dict:
                print('in table')
                group_dict[suf].append(key)
            else:
                group_dict[suf] = [key]
        else:
            group_dict[key] = [key]


    #plot all calculated effective variables
    for var,vars in group_dict.items():
        plt.clf()
        print(f'var = {var}')
        for key in vars: #iterate over the output gating parameters
            print(f'key = {key}')
            if 'tcomp' in key:
                continue #don't plot the computation time tcomp as these plots are not important
            table = pkldict['tables'][key] #gating parameter matrix, dims:(a,f,A,Q,Cm,fs)
            Q_array = Qrange if len(Qrange) in table.shape else Q_ext #if the length of Qrange is the same as on of the dimensions in the table, use this array, otherwise use the extended Q array
            table2 = eval(f"table[{(',').join(ind_list)}]") #index the table according to the given indexing list: 1 value for each parameter except all Q and a chosen variable
            if all_list:
                ind_colon = np.where(np.array(ind_list)==":")[0] #determine which variables keep all values (Q and a chosen variable)
                ind = ind_colon[ind_colon !=3][0] #except from Q (which has index 3), check the index of the variable that has all its values
                if ind > 3: #if the variable is further in the array than Q:
                    table2 = np.swapaxes(table2,0,1) #switch the variable with the Q-axis so Q is plotted on the x-axis and the other one on the y-axis as desired
                table2 = np.reshape(table2,(-1,len(Q_array))) #reshape the table to dim(variable) x dim(Q) (this is partly done because of fs-dimension)
            plt_LUTeff(key,table2,Q_array,factor,unit,var_refs=var_refs,reduced_yrange=reduced_yrange,reduced_xrange=reduced_xrange, label = 1) #plotting
            #plt.title(title)
        plt.legend()
        plt.savefig(f'figs/{foldername}/{key}.png')

    #plot the effective capacity
    table_V = eval(f"pkldict['tables']['V'][{(',').join(ind_list)}]") #index the table according to the calculated indexing list
    if all_list:
        if ind > 3: #if the variable is further in the array than Q:
            table_V = np.swapaxes(table_V,0,1) #switch the variable with the Q-axis so Q is plotted on the x-axis and the other one on the y-axis as desired
        table_V = np.reshape(table_V,(-1,len(Qrange))) #reshape the table to dim(variable) x dim(Q)
    #Q_table = np.tile(Qrange,np.prod(table_V.shape)//len(Qrange)).reshape(table_V.shape) #copy the Q array to a matrix with the same dimensions as V -> not needed, use broadcasting
    table_C = Qrange / table_V #Q = C x V <=> C = Q / V: C/ m2 / mV = 1 / 1e-3 * C/m2/V = 1e3 F/m2 = 1e5 uF/m2
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
        plt.clf()
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


def plot_astim(csv_file, separate=False, debug = False, folder = r'C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 7\\'):
    """to plot all the plotting variables (including gating parameters) in a certain compartment of a time simulation
        :csv_file: file containing the time variable and all variables during the timelapse
        :separate: plot all the variables also on separate plots
        :folder: directory where the plot is saved """
    
    with open(csv_file, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            titles = row[1:] #first string in header is empty (row index)
            break
    sim_csv = np.genfromtxt(csv_file, delimiter=',',skip_header=1)
    sim_csv = sim_csv.swapaxes(1,0)[1:] #first column (now row after swapping) contains row index
    titles = [e.split(r"'")[1]+"_"+e.split(r"'")[3] if ('(' in e and ')' in e) else e for e in titles] #the gating parameters are stored as (x, mech), convert it to: x_mech
    factors = [1e3, 1, 1e5, 1] #factors for time, stimstate, membrane  charge, and voltage(which is already in mV)
    labels = ['t [ms]', 'stimstate', 'Q [nC/cm2]', 'V [mV]'] 
    factors_add = [1 if (e.startswith('i') or e.startswith('I')) else 1e2 if 'Cm' in e else 1 if ('_' in e) else 1 for e in titles[4:]] #the factors are added for the gating parameters, membrane capacitance and currents
    #factors_add = [1 if ('_' in e) else 1e2 if 'Cm' in e else 1 if ('i' in e or 'I' in e) else 0 for e in titles[4:]] #the factors are added for the gating parameters, membrane capacitance and currents
    labels_add = [f'{e} [mA/m2]' if (e.startswith('i') or e.startswith('I')) else 'Cm [uF/cm2]' if 'Cm' in e else f'{e}' if ('_' in e) else '' for e in titles[4:]] #labels are added for the gating parameters, membrane capacitance and currents
    #labels_add = [f'{e}' if ('_' in e) else 'Cm [uF/cm2]' if 'Cm' in e else f'{e} [mA/m2]' if ('i' in e or 'I' in e) else 0 for e in titles[4:]] #labels are added for the gating parameters, membrane capacitance and currents
    factors += factors_add
    labels += labels_add
    print(f"sim step: {sim_csv[0][2]*1e6} us, end sim: {sim_csv[0][-1]*1e3:.3f} ms")
    if separate: #to plot every variable separately
        for i in range(0,len(sim_csv)-1): #len(sim_csv)
            plt.plot((sim_csv[0]*factors[0])[:], (sim_csv[i+1]*factors[i+1])[:]) #slicing until x for debugging
            plt.title(titles[i+1])
            plt.xlabel(labels[0])
            plt.ylabel(labels[i+1])
            plt.show()
    plot_QmVm = 1
    if plot_QmVm: #plot Qm/Vm which should be Cm (with or without a delay)
        plt.plot((sim_csv[0]*factors[0])[:], (sim_csv[2]*factors[2]/sim_csv[7]*factors[7])[:])
        plt.show()
    nrows = int(np.ceil((len(sim_csv)-1)/2)) #number of columns = 2, number of rows depends on the number of variables
    fig1, axs = plt.subplots(nrows, 2, figsize=(15,7))
    for i in range(0,len(sim_csv)-1): #len(sim_csv)
        if 'nC/cm2' in labels[i+1]:
            Qrow = (sim_csv[i+1]*factors[i+1])
        if 'mV' in labels[i+1]:
            Vrow = (sim_csv[i+1]*factors[i+1])
        if 'uF/cm2' in labels[i+1]:
            Crow = (sim_csv[i+1]*factors[i+1])
        trow = sim_csv[0]
        axs[i//2,i%2].plot((sim_csv[0]*factors[0])[:], (sim_csv[i+1]*factors[i+1])[:])
        axs[i//2,i%2].set_ylabel(labels[i+1],fontsize=8,rotation=0,labelpad=30)
        if i//2 == nrows-1:
            axs[i//2,i%2].set_xlabel(labels[0]) #only set x label at the bottom
            axs[i//2,i%2].set_xticks((sim_csv[0]*factors[0])[::len(sim_csv[0])//9]) #so time is plotted correctly on x-axis #only set x ticks at the bottom
        else:
            axs[i//2,i%2].set_xticks([])
        #axs[i//2,i%2].set_title(titles[i+1])
    if debug:
        arg = np.argmin(abs(Vrow))
        print(arg/len(Vrow))
        #print((Qrow)[arg-2:arg+2])
        #print((Vrow)[arg-2:arg+2])
        #print(Crow[arg-2:arg+2])
        # print((Qrow)[:5])
        # print((Vrow)[:5])
        # print(Crow[:5])
        # print(trow[:5])
    axs[nrows-1,1].set_xlabel(labels[0]) #do this again in case there is no curve plotted on the last subplot
    axs[nrows-1,1].set_xticks((sim_csv[0]*factors[0])[::len(sim_csv[0])//9]) #so time is plotted correctly on x-axis #in case of odd number of plots as explained in line above
    plt.subplots_adjust(hspace=.0)
    fig1.suptitle("Stim")
    fig1.tight_layout()
    path_pieces = csv_file.split('\\') #path is splitted into its directories
    directory = folder+path_pieces[-2]+'_ext\\'
    if not os.path.exists(directory):
        os.mkdir(directory)  

    if debug:
        #quit()
        plt.show()
    if not debug:
        plt.savefig(directory+path_pieces[-1].replace('.csv',f'.jpg')) #save the image
        print(f"Image saved at: {directory+path_pieces[-1].replace('csv','jpg')}")


def plot_astim2(csv_file, separate=0, debug=0, variables=None, folder=None, save_fig = 0, save_sv=0, tauinf=0, compare=0):#, folder = r'C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 7\\'):
    """to plot all the plotting variables (including gating parameters) in a certain compartment of a time simulation
        :csv_file: file containing the time variable and all variables during the timelapse
        :separate: plot all the variables also on separate plots
        :variables: simulated variables that need to be plotted
        :folder: directory where the plot is saved """
    
    sim_csv, titles, factors, labels = read_csv(csv_file)

    if save_sv:
        #to save a certain state vector
        pre_offset = np.argmin(abs(sim_csv[0] - (1)*1e-3))
        post_offset = pre_offset+1
        pre = [sim_csv[i][pre_offset] for i in range(len(sim_csv))]
        post = [sim_csv[i][post_offset] for i in range(len(sim_csv))]
        print(pre,'\n', post,'\n', np.array(pre)-np.array(post))
        quit()

    path_pieces = csv_file.split('\\') #path is splitted into its directories
    if folder:
        directory = folder+path_pieces[-2]+'_ext\\'
        if not os.path.exists(directory):
            os.mkdir(directory)  

    vars = [e.split(' ')[0] for e in labels]
    variables = variables if variables else copy.copy(vars)
    variables.remove('t') if 't' in variables else None

    plot_dict ,Qrow, Vrow, loc = create_plot_dict(titles,factors,labels,vars,variables,sim_csv)
    #plot_dict['dQ/dt'] = {'label': 'dQ/dt [mA/m2]', 'array': np.diff(plot_dict['Q']['array'][1:])/np.diff(sim_csv[0][1:]), 'title': 'dQ/dt', 'factor': 1e3, 'loc' : loc} #remove the first value as these are the same as the second resulting in 0/0
    #plot_dict['dQ/dt']['array'] = np.append(np.append([0],(plot_dict['dQ/dt']['array'])),[0]) #np.append(plot_dict['dQ/dt']['array'][1:],plot_dict['dQ/dt']['array'][-2:])    

    "only do this in case of backward Euler which has a delay of 1 sample"
    # if 'Cm' in plot_dict:
    #     Crow = Qrow / np.concatenate((Vrow[1:], [Vrow[-1]])) / 1e2 # uF/cm2 -> F/m2
    #     plot_dict['Cm']['array'] = Crow

    if tauinf:
        plot_dict['V']['array'] = np.linspace(-150,50,400) #in case tau and inf need to be plotted in function of V instead of t
        time =  np.linspace(-150,50,400) #(sim_csv[0]*factors[0])[:] #options: in function of V or t
        mtau_SKv3_1 = 0.2*20.000/(1+np.exp(((plot_dict['V']['array'] -(-46.560))/(-44.140))))
        minf_SKv3_1 = 1/(1+np.exp(((plot_dict['V']['array'] -(18.700))/(-9.700))))

        malpha_NaTs2_t = (0.182 * (plot_dict['V']['array']- -32))/(1-(np.exp(-(plot_dict['V']['array']- -32)/6)))
        mbeta_NaTs2_t = (0.124 * (-plot_dict['V']['array'] -32))/(1-(np.exp(-(-plot_dict['V']['array'] -32)/6)))
        halpha_NaTs2_t = (-0.015 * (plot_dict['V']['array']- -60))/(1-(np.exp((plot_dict['V']['array']- -60)/6)))
        hbeta_NaTs2_t = (-0.015 * (-plot_dict['V']['array'] -60))/(1-(np.exp((-plot_dict['V']['array'] -60)/6)))
        qt = 2.3**((37-21)/10)

        mtau_NaTs2_t = (1/(malpha_NaTs2_t + mbeta_NaTs2_t))/qt
        minf_NaTs2_t = malpha_NaTs2_t/(malpha_NaTs2_t + mbeta_NaTs2_t)
        htau_NaTs2_t = (1/(halpha_NaTs2_t + hbeta_NaTs2_t))/qt
        hinf_NaTs2_t = halpha_NaTs2_t/(halpha_NaTs2_t + hbeta_NaTs2_t)
        plt.plot(time, mtau_SKv3_1*1e2,label='m_tau_SKv3_1')
        plt.plot(time, minf_SKv3_1*1e2,label='m_inf_SKv3_1')
        plt.plot(time, mtau_NaTs2_t*1e2,label='m_tau_NaTs2_t')
        plt.plot(time, minf_NaTs2_t*1e2,label='m_inf_NaTs2_t')
        plt.plot(time, htau_NaTs2_t*1e2,label='h_tau_NaTs2_t')
        plt.plot(time, hinf_NaTs2_t*1e2,label='h_inf_NaTs2_t')
        plt.legend()
        plt.grid()
        plt.show()
        quit()

    if compare: #plot a selection of variables (scaled) on top of each other (on a single plot)
        for var in plot_dict.keys():#variables:
            #print(var)
            if var == 'Q' or var == 'V':
                continue
                plt.plot((sim_csv[0]*factors[0])[:], (plot_dict[var]['array']*plot_dict[var]['factor'])[:],label=var) #slicing until x for debugging
            elif var == 'i_net':# or var == 'i_SKv31' or var == 'i_NaTs2t':
                pass
                plt.plot((sim_csv[0]*factors[0])[:], (plot_dict[var]['array']*plot_dict[var]['factor'])[:],label=var) #slicing until x for debugging
            elif 'SK' in var or 'Na' in var:
                continue
                plt.plot((sim_csv[0]*factors[0])[:], (10**2.3)*abs(plot_dict[var]['array']*plot_dict[var]['factor'])[:],label=var) #slicing until x for debugging
        print(-np.diff(plot_dict['Q']['array'][1:])/np.diff(sim_csv[0][1:]))
        print(plot_dict['Q']['factor'])
        plot_dict['dQ/dt'] = {'label': 'dQ/dt [mA/m2]', 'array': np.diff(plot_dict['Q']['array'][1:])/np.diff(sim_csv[0][1:]), 'title': 'dQ/dt', 'factor': 1e3, 'loc' : loc} #remove the first value as these are the same as the second resulting in 0/0
        plot_dict['dQ/dt']['array'] = np.append(np.append([0],(plot_dict['dQ/dt']['array'])),[0]) #np.append(plot_dict['dQ/dt']['array'][1:],plot_dict['dQ/dt']['array'][-2:])    
        plt.plot((sim_csv[0]*factors[0])[:],-plot_dict['dQ/dt']['array']*plot_dict['Q']['factor']*1e3,label = 'dQ/dt')
        plt.legend()
        plt.grid()
        plt.show()        
        quit()
    if debug:
        #plt.plot(np.diff(sim_csv[0])*1e3); plt.show()
        #arg = np.argmin(abs(Vrow))
        #print(arg/len(Vrow))
        #print((Qrow)[arg-2:arg+2]); print((Vrow)[arg-2:arg+2]); print(Crow[arg-2:arg+2])
        # print((Qrow)[:5]); print((Vrow)[:5]); print(Crow[:5]); print(trow[:5])
        quit()
    if separate: #to plot every variable separately
        for var in plot_dict.keys():#variables:
            plt.clf()
            plt.plot((sim_csv[0]*factors[0])[:], (plot_dict[var]['array']*plot_dict[var]['factor'])[:],label=' section ') #slicing until x for debugging
            #plt.axvline(x=108,color='r') #plot a vetical line when stimulation ends
            #plt.title(''plot_dict[var]['title']'')
            plt.xlabel(labels[0])
            plt.ylabel(plot_dict[var]['label'])
            plt.legend()
            plt.grid()
            if folder:
                print(directory)
                var_pure = var.replace('\\','').replace('/','')
                plt.savefig(directory+r'\\'+f'{var_pure}.jpg')
            else:
                plt.show()
        quit()    
    
    nrows = int(np.ceil((loc+1)/2)) #number of columns = 2, number of rows depends on the number of variables # no -1 because variables doesn't contain the time (x-axis) as in plot_astim(1)
    fig1, axs = plt.subplots(nrows, 2, figsize=(15,7))
    print(f'figure created with dimensions: {nrows}, 2')
    for i,var in enumerate(plot_dict.keys()): #len(sim_csv)
        var_dict = plot_dict[var]
        trow = sim_csv[0]
        axs[var_dict['loc']//2,var_dict['loc']%2].plot((sim_csv[0]*factors[0])[:], (var_dict['array']*var_dict['factor'])[:],label = var_dict['label'])
        axs[var_dict['loc']//2,var_dict['loc']%2].set_ylabel(var_dict['label'],fontsize=8,rotation=0,labelpad=30)
        if var_dict['loc']//2 == nrows-1:
            axs[var_dict['loc']//2,var_dict['loc']%2].set_xlabel(labels[0]) #only set x label at the bottom
            #axs[var_dict['loc']//2,var_dict['loc']%2].set_xticks((sim_csv[0]*factors[0])[::len(sim_csv[0])//9]) #so time is plotted correctly on x-axis #only set x ticks at the bottom
        elif loc > 3:
            axs[var_dict['loc']//2,var_dict['loc']%2].set_xticks([])
        #axs[i//2,i%2].set_title(titles[i+1])      
    axs[nrows-1,1].set_xlabel(labels[0]) #do this again in case there is no curve plotted on the last subplot
    #axs[nrows-1,1].set_xticks((sim_csv[0]*factors[0])[::len(sim_csv[0])//9]) #so time is plotted correctly on x-axis #in case of odd number of plots as explained in line above
    plt.subplots_adjust(hspace=.0)
    fig1.suptitle("Stim")
    fig1.tight_layout()

    for i in range(nrows):
        for j in range(2):
            axs[i,j].legend(loc='upper left', fontsize=5) #left when doing stimulation, right when debugging
            axs[i,j].grid()
    if save_fig:
        plt.savefig(directory+path_pieces[-1].replace('.csv',f'.jpg')) #save the image
        print(f"Image saved at: {directory+path_pieces[-1].replace('csv','jpg')}")
    else:
        plt.show()


def plot_astim_sections(csv_files,variables=None, debug=0):
    plot_dicts = {}
    for csv_file in csv_files:
        sim_csv, titles, factors, labels = read_csv(csv_file)
        
        vars = [e.split(' ')[0] for e in labels]
        variables = variables if variables else copy.copy(vars)
        variables.remove('t') if 't' in variables else None
        plot_dict ,Qrow, Vrow, loc = create_plot_dict(titles,factors,labels,vars,variables,sim_csv)
        plot_dicts[csv_file.split('_')[-1]] = plot_dict
    for e,f in plot_dicts.items():
        var_plot = 'V'
        x, y = 0, -1
        #dict_keys(['stimstate', 'Q', 'V', 'Cm', 'iax', 'i_net', 'i_pas'])
        # if 'node' in e:
        #     var_dict = f[var_plot]
        #     plt.plot((sim_csv[0]*factors[0])[x:y], (var_dict['array']*var_dict['factor'])[x:y],label = e)
        #     a = var_dict['array']*var_dict['factor']
        # elif e.startswith('myelin'):
        #     var_dict = f[var_plot]
        #     plt.plot((sim_csv[0]*factors[0])[x:y], (var_dict['array']*var_dict['factor'])[x:y],label = e)
        #     b = var_dict['array']*var_dict['factor']
        # else: continue
        var_dict = f[var_plot]
        color = 'yellow' if 'soma' in e else 'purple' if 'axon' in e else '#FF0000' if 'node' in e else '#FFA500' if 'unmyelin' in e else '#000000' if 'myelin' in e else '#0000FF' if 'apical' in e else '#00FF00' if 'basal' in e else None
        plt.plot((sim_csv[0]*factors[0])[:], (var_dict['array']*var_dict['factor'])[:],label = e.split('0.csv')[0], color=color)
    if var_plot:
        pass
        #plt.title(var_plot)
    if debug:
        print(a[0],b[0])
        print(np.mean(a/b),np.max(a/b),np.min(a/b))
    fs = 36
    plt.xlabel('time [ms]',fontsize=fs)
    plt.ylabel('membrane voltage [mV]',fontsize=fs)
    plt.xticks(fontsize = fs)
    plt.yticks(fontsize = fs)
    plt.legend(fontsize = fs-4,loc='lower right')
    plt.grid()
    plt.show()


def plot_titration_curves(pklfile):
    """to plot titration curve in function of the duty cycle"""
    pkldict = read_pickle(pklfile)
    PRF_var = 0
    freq_var = 0
    rad_var = 0
    fs = 20

    plt.rc('xtick',labelsize=fs)
    plt.rc('ytick',labelsize=fs)

    #pprint.pprint(pkldict)
    #print('cell freq\tradius\t fs\tDC\tPRF\tamp')
    #pprint.pprint(pkldict)
    cell_nr = 7
    freq = 500000.0 #500000.0 #100000.0
    freq1 = 100000.0
    freq2 = 500000.0
    freq3 = 1000000.0
    radius = 3.2e-08 #3.2e-08 #1.6e-08
    radius1 = 1.6e-08
    radius2 = 3.2e-08
    radius3 = 6.4e-08
    coverage = 0.75
    DC = all
    PRF = 1000.0
    PRF1 = 10.0
    PRF2 = 100.0 #10.0 #50.0
    PRF3 = 1000.0

    #x = np.arange(0.05,1.05,0.05)#[1:]
    #x = np.arange(0.1,0.67,0.01) #np.arange(0.1,1.01,0.01)
    x = np.array(range(15,101))/100
    x2 = np.array(range(15,101,5))/100
    DC = np.array(range(15,101))/100
    DC1, DC2, DC3 = DC, DC, DC
    #print(x2)
    #print(pkldict[7].keys())
    if PRF_var:
        TIT1 = np.array([pkldict[cell_nr][freq][radius][coverage][e][PRF1] for e in DC1])
        TIT2 = np.array([pkldict[cell_nr][freq][radius][coverage][e][PRF2] for e in DC2])
        TIT3 = np.array([pkldict[cell_nr][freq][radius][coverage][e][PRF3] for e in DC3])
        #filtering of unexcited values
        DC1, DC2, DC3 = DC1[TIT1>0], DC2[TIT2>0], DC3[TIT3>0]
        TIT1, TIT2, TIT3 = TIT1[TIT1>0], TIT2[TIT2>0], TIT3[TIT3>0]        
        plt.plot(DC1,TIT1*1e-3,label=str(PRF1)+' Hz',color='cyan')
        plt.plot(DC2,TIT2*1e-3,label=str(PRF2)+' Hz',color='deepskyblue')
        plt.plot(DC3,TIT3*1e-3,label=str(PRF3)+' Hz',color='royalblue')
    if freq_var:
        TIT1 = np.array([pkldict[cell_nr][freq1][radius][coverage][e][PRF] for e in DC1])
        TIT2 = np.array([pkldict[cell_nr][freq2][radius][coverage][e][PRF] for e in DC2])
        TIT3 = np.array([pkldict[cell_nr][freq3][radius][coverage][e][PRF] for e in DC3])
        #filtering of unexcited values
        DC1, DC2, DC3 = DC1[TIT1>0], DC2[TIT2>0], DC3[TIT3>0]
        TIT1, TIT2, TIT3 = TIT1[TIT1>0], TIT2[TIT2>0], TIT3[TIT3>0]
        plt.semilogy(DC1,TIT1*1e-3,label=str(freq1*1e-3)+' kHz',color='cyan')
        plt.semilogy(DC2,TIT2*1e-3,label=str(freq2*1e-3)+' kHz',color='deepskyblue')
        plt.semilogy(DC3,TIT3*1e-3,label=str(freq3*1e-3)+' kHz',color='royalblue')
    if rad_var:
        TIT1 = np.array([pkldict[cell_nr][freq][radius1][coverage][e][PRF] for e in DC1])
        TIT2 = np.array([pkldict[cell_nr][freq][radius2][coverage][e][PRF] for e in DC2])
        TIT3 = np.array([pkldict[cell_nr][freq][radius3][coverage][e][PRF] for e in DC3])
        #filtering of unexcited values
        DC1, DC2, DC3 = DC1[TIT1>0], DC2[TIT2>0], DC3[TIT3>0]
        TIT1, TIT2, TIT3 = TIT1[TIT1>0], TIT2[TIT2>0], TIT3[TIT3>0]
        plt.semilogy(DC1,TIT1*1e-3,label=str(radius1*1e9)+' nm',color='cyan')
        plt.semilogy(DC2,TIT2*1e-3,label=str(radius2*1e9)+' nm',color='deepskyblue')
        plt.semilogy(DC3,TIT3*1e-3,label=str(radius3*1e9)+' nm',color='royalblue')
    else:
        TIT_info = np.array([pkldict[cell_nr][freq][radius][coverage][e][PRF] for e in DC])
        TIT = np.array([pkldict[cell_nr][freq][radius][coverage][e][PRF][0]*1e-3 for e in DC])
        #filtering of unexcited values
        DC = DC[TIT>0]
        TIT = TIT[TIT>0]
        labels = []
        for i,e in enumerate(TIT_info):
            #print(e)
            color = 'yellow' if 'soma' in e[1] else 'purple' if 'axon' in e[1] else '#FF0000' if 'node' in e[1] else '#FFA500' if 'unmyelin' in e[1] else '#000000' if 'myelin' in e[1] else '#0000FF' if 'apical' in e[1] else '#00FF00' if 'basal' in e[1] else None
            label = 'soma' if 'soma' in e[1] else 'axon' if 'axon' in e[1] else 'node' if 'node' in e[1] else 'unmyelin' if 'unmyelin' in e[1] else 'myelin' if 'myelin' in e[1] else 'apical' if 'apical' in e[1] else 'basal' if 'basal' in e[1] else None
            if label in labels:
                label = None
            else:
                labels.append(label)
            print(f'{cell_nr},{int(radius*1e9)},{int(freq*1e-3)},{int(PRF)},{int(DC[i]*1e2)},{TIT[i]+0.0390625}')
            plt.scatter(DC[i],TIT[i],color=color,label = label)
        plt.yscale('log')
        #plt.semilogy(DC,TIT*1e-3,label=str(radius1*1e9)+' nm',color='cyan')
    plt.xlabel('Duty cycle [%]',fontsize=fs)
    plt.ylabel('Amplitude [kPa]',fontsize=fs)
    plt.xticks(fontsize = fs)
    plt.yticks(fontsize = fs)
    plt.legend(fontsize = fs)
    plt.show()


def plot_titration_curves_cells(pklfile):
    """to plot titration curve in function of the duty cycle"""
    pkldict = read_pickle(pklfile)
    #pprint.pprint(pkldict)
    #print('cell freq\tradius\t fs\tDC\tPRF\tamp')
    #pprint.pprint(pkldict)
    cell_nr = 7
    freq = 500000.0 #500000.0 #100000.0
    radius = 3.2e-08 #3.2e-08 #1.6e-08
    coverage = 0.75
    DC = all
    PRF1 = 100.0
    PRF2 = 10.0 #10.0 #50.0
    PRF3 = 1000.0
    labels = ["L1 NGC-DA", "L2/3 PC", "L4 LBC", "L5 TTPC", "L6 TPC"]
    
    for freq in [500000.0, 1000000.0]:
        for radius in [3.2e-8]:#, 6.4e-8]:
            for PRF in [10.0, 100.0, 1000.0]:
                for cell_group in range(5):
                    bp = np.empty((5,2)) #we create a 5x2 matrix which will plot for the 5 cell groups the value of 0.5 and 1.0 for DC
                    for i,cell_nr in enumerate(np.array([e for e in range(5)]) + cell_group*5 + 1): #iterate over the 5 cells in each group
                        DC = [0.5, 1.0]
                        TIT = [pkldict[cell_nr][freq][radius][coverage][e][PRF] for e in DC]
                        TIT = [e if type(e) == float else e[0] for e in TIT]
                        TIT = np.array(TIT)
                        bp[i] = TIT*1e-3
                        #plt.plot(DC,TIT*1e-3,label=cell_nr)
                    plt.plot(DC,np.median(bp,axis=0),label=labels[cell_group])
                    plt.fill_between(DC,np.quantile(bp,0.45,axis=0),np.quantile(bp,0.55,axis=0),alpha=0.3)
                    #plt.boxplot(bp,positions=range(2),labels=DC)
                fs = 18
                plt.xlabel('Duty cycle [%]',fontsize=fs)
                plt.ylabel('Amplitude [kPa]',fontsize=fs)
                plt.xticks(fontsize = fs)
                plt.yticks(fontsize = fs)
                plt.legend(fontsize = fs)
                plt.title(f'freq = {freq}, radius = {radius}, PRF = {PRF}')
                plt.show()


"""------------------------------------------------------------------------------------ANALYSIS------------------------------------------------------------------------------------"""
def analyze_over_sections(pkl_file,pkl_file2=None):
    """takes the whole dictionary (pkl file) of a simulation and analyzes different metrics accross all sections and also for a specific section"""

    pkldict = read_pickle(pkl_file)
    if pkl_file2:
        pkldict2 = read_pickle(pkl_file2)
    section = 'soma0' #considered section that will be plotted
    t_mini = 1e10
    sections = [] #all sections
    t_first_zero_crossings = [] #all time points where first zero-crossing happens for every section
    V_first_zero_crossings = [] #voltage value before the first zero-crossing (can also be the one after because two sections aren't going to cross each other during the AP, right? -> THEY CAN)
    V0, V0_2 = [],[]

    for e,f in pkldict.items(): #we iterate over all sections and their dictionary
        crossings = ((np.array(f['Vm'])[:-1] * np.array(f['Vm'])[1:]) <= 0)*1 #V/voltage zero-crossings of section e (0 if no crossing, 1 if crossing)
        first_crossing = np.argmax(crossings)+1 if np.argmax(crossings) != 0 else len(crossings) #index value of first zero-crossing
        t_zero_crossing = f['t'][first_crossing] #timestamp where first zero-crossing happens
        V_zero_crossing = f['Vm'][first_crossing-1] #voltage value before the zero-crossing
        #print(e,np.min(abs(f['Vm'])))
        sections.append(e) #append every section to the sections list
        if section in e: #if we are at the specified section: 
            V_section, t_section = np.array(f['Vm']), np.array(f['t']) #store the t and V array of this section
            t_stim = t_section[np.array(f['stimstate']>0)] #store the stimulation time (this is the time array when the stimulation is on)
        t_first_zero_crossings.append(t_zero_crossing) #all the zero-crossings timestamps for every section are stored
        V_first_zero_crossings.append(V_zero_crossing) #all the zero-crossings voltages for every section are stored
        #V0.append(f['Vm'][0]) #append Vm0 for every section (too analyze what the init value is)
        if pkl_file2:
            pass
            #V0_2.append(pkldict2[e]['Vm'][0])
        if t_zero_crossing < t_mini: #here we determine what the smallest first zero-crossing is
            t_mini = t_zero_crossing
        #color scheme of aberra cells
        color = 'yellow' if 'soma' in e else 'purple' if 'axon' in e else '#FF0000' if 'node' in e else '#FFA500' if 'unmyelin' in e else '#000000' if 'myelin' in e else '#0000FF' if 'apical' in e else '#00FF00' if 'basal' in e else None
        #label every type of section by using the label for the section0
        #continue
        if (e[-1] == '0' and e[-2].isalpha()):
            if pkl_file2:
                # print(f['t']*1e3)
                # print(pkldict2[e]['t']*1e3)
                # print(np.append(f['t'],pkldict2[e]['t']*1e3))
                # quit()
                if e == 'axon':
                    e = 'axon initial segment'
                plt.plot(np.append(f['t']*1e3,(pkldict2[e]['t']+np.array(f['t'])[-1])*1e3),np.append(f['Vm'],pkldict2[e]['Vm']),label=e[:-1], color=color)
            else:
                #print(e)
                if e[:-1] == 'axon':
                    e = 'AIS0'
                plt.plot(f['t']*1e3,f['Vm'],label=e[:-1], color=color)
        else:
            if pkl_file2:
                plt.plot(np.append(f['t']*1e3,(pkldict2[e]['t']+np.array(f['t'])[-1])*1e3),np.append(f['Vm'],pkldict2[e]['Vm']), color=color)
            else:
                plt.plot(f['t']*1e3,f['Vm'], color=color)

    sections_sorted = [section for _, _, section in sorted(zip(t_first_zero_crossings, abs(np.array(V_first_zero_crossings)) ,sections))] #first we look at the time that took until the first zero-crossing happened
    zero_crossing_vals_sorted = [val for _, _, val in sorted(zip(t_first_zero_crossings, abs(np.array(V_first_zero_crossings)), V_first_zero_crossings))] #then we look at the section that is closest to zero before crossing it
    act_section = sections_sorted[0]
    # print(sections)
    # print(zero_crossings)
    fs = 30
    plt.xlabel('time [ms]',fontsize=fs)
    plt.ylabel('membrane voltage [mV]',fontsize=fs)
    plt.xticks(fontsize = fs)
    plt.yticks(fontsize = fs)
    plt.legend(fontsize = fs)
    plt.show()
    #plt.plot(t_section*1e3,V_section)
    # plt.show()

    #metrics
    #nAPs = len(np.where(np.diff(np.sign(V_soma)))[0])/2 #number of action potentials by looking at the number of zero-crossings/2
    pos_crossings = (np.array(V_section)[:-1] <= 0)*(np.array(V_section)[1:] > 0) #crossing where the value goes from a negative to a positive one
    pos_crossings = np.append([0],pos_crossings) #add 0 to comply with the time array length

    spiking_times = t_section*pos_crossings
    spiking_times = spiking_times[spiking_times!=0] #times where the V value is positive after going through a zero-crossing

    #Vmax_index = argrelextrema(np.array(V_section), np.greater)[0] #looks at relative maxima
    #Vmax_t = t_section[Vmax_index] #time at relative maxima
    #Vmax_V = V_section[Vmax_index] #voltage at relative maxima
    #spiking_times = Vmax_t[Vmax_V>0] #we assume that there is a spike if V is at a local maxima above zero -> very wrong assumption
    spiking_and_bounds = np.append(np.append(t_stim[0],spiking_times),t_stim[-1])
    ISI = np.diff(spiking_times)*1e3
    ISI_with_bounds = np.diff(spiking_and_bounds)*1e3
    for e in spiking_times:
        pass
        #plt.axvline(x=e*1e3,c='red')
    nAPs = len(spiking_times)

    print(f'(number of sections: {len(sections)})')
    print(f'activation section/site:{act_section} @ {t_mini*1e3} ms') #section where the first zero-crossing happens so first AP
    for e,f,g in zip(sections_sorted[:10],sorted(t_first_zero_crossings)[:10],zero_crossing_vals_sorted[:10]):
        print(f'{e}: \t{f},\t{g}')
    #print(f'V0: {V0[:5]}, V0_2: {V0_2[:5]}')
    print(f"number of APs: {nAPs}")
    print(f"ISI: {ISI}")
    firing_rate = 1/ISI*1e3
    #plt.plot(np.cumsum(ISI),firing_rate)
    print(f"ISI with bounds: {ISI_with_bounds}")
    #print(f"spiking times: {spiking_times}")
    #print([x for _, x in sorted(zip(zero_crossings, sections))]) #sections ordered by first stimulation moment
    #plt.show()
   

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