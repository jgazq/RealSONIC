# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-09-27 14:28:52
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2023-03-29 18:39:34

''' Run simulations of a multi-compartmental realistic Aberra/BBP neuron model with a specific point-neuron mechanism
    upon ultrasound stimulation at one node. '''

from PySONIC.core import Batch, NeuronalBilayerSonophore
from PySONIC.utils import logger, gaussian, gaussian3D, loadData
from PySONIC.parsers import AStimParser
from MorphoSONIC.core import GaussianAcousticSource, UniformAcousticSource, Gaussian3DAcousticSource
from MorphoSONIC.parsers import AStimRealisticNeuronParser

import pickle
import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import os
from neuron import h
import sys


def simulate_realnrn(args,amp,stimulation,x0,sigma, fiber, iter):

    output = []
    args['amp'] = [amp]

    queue = [item[:2] for item in NeuronalBilayerSonophore.simQueue(
        *AStimParser.parseSimInputs(args), outputdir=args['outputdir'])]
    if args['save']:
        queue = [(item[0][:2], item[1]) for item in queue]
   
    if args['save']:
        if stimulation == "gaussian":
            simqueue = [(
                [GaussianAcousticSource(
                    x0, sigma, item[0][0].f, item[0][0].A), *item[0][1:]], #[1:]: the AcousticDrive is replaced by and AcousticSource
                item[1]                                                    #both Drive and Source inherit from StimObject
            ) for item in queue]
        elif stimulation == "uniform":
            simqueue = [(
                [UniformAcousticSource(
                    item[0][0].f, item[0][0].A), *item[0][1:]],
                item[1] #here you can add the arguments of dt and atol
            ) for item in queue]   
        elif stimulation == "gaussian3D":
            simqueue = [(
                [Gaussian3DAcousticSource(
                    x0, x0, x0, sigma, sigma, sigma, item[0][0].f, item[0][0].A), *item[0][1:]], #currently the same x0 and sigma is given in all 3 directions but this should be variable
                item[1]
            ) for item in queue]                                                                
        func = fiber.simAndSave
    else:
        if stimulation == "gaussian":
            simqueue = [
                [GaussianAcousticSource(x0, sigma, item[0].f, item[0].A), *item[1:]]
                for item in queue]
        elif stimulation == "uniform":
            simqueue = [
                [UniformAcousticSource(item[0].f, item[0].A), *item[1:]]
                for item in queue]   
        elif stimulation == "gaussian3D":     
            simqueue = [
                [Gaussian3DAcousticSource(x0, x0, x0, sigma, sigma, sigma, item[0].f, item[0].A), *item[1:]]
                for item in queue]                                                        
        func = fiber.simulate
    batch = Batch(func, simqueue)
    output += batch(loglevel=args['loglevel']) #BREAKPOINT
    #print(args['section'])
    refsec = fiber.sections[args['section'][0][:-1]][args['section'][0]] #fiber.refsection
    amplitude = gaussian(refsec.x_xtra,x0,sigma,queue[0][0].A) if stimulation == 'gaussian' else gaussian3D(refsec.x_xtra, refsec.y_xtra, refsec.y_xtra,x0,x0,x0,sigma,sigma,sigma,queue[0][0].A) if stimulation == 'gaussian3D' else queue[0][0].A
    #print(f'{"-"*50}\nrefsection:\nlocation:\t({refsec.x_xtra}, {refsec.y_xtra}, {refsec.z_xtra}) um\namplitude = {amplitude*1e-3}kPa')
    args['ref_loc'] = (refsec.x_xtra, refsec.y_xtra, refsec.z_xtra)

    tosave = output[0][0].data
    outdir = r"C:\Users\jgazquez\RealSONIC\pickledump" + f"\{args['fs'][0]*100}%_{args['radius'][0]*1e9}nm_{args['freq'][0]*1e-3}kHz"  + f"_{args['tstim'][0]*1e3}ms_{args['toffset'][0]*1e3}ms_{args['PRF'][0]}Hz_{args['DC'][0]}DC\\"
    outname = r"dump_" + f"{args['fs'][0]*100}%_{args['radius'][0]*1e9}nm_{args['freq'][0]*1e-3}kHz"  + f"_{args['tstim'][0]*1e3}ms_{args['toffset'][0]*1e3}ms_{args['PRF'][0]}Hz_{args['DC'][0]}DC_{args['amp'][0]*1e-3}kPa" + ".pkl"
    #outfile = outdir + outname
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        os.mkdir(outdir+'pkl\\')
        os.mkdir(outdir+'csv\\')

    with open(outdir+'pkl\\'+outname, 'wb') as fh:
        pickle.dump(tosave, fh)
    data = tosave['soma0']
    data.to_csv((outdir+'csv\\'+outname).replace('pkl','csv'))
    return tosave

def thresh_excited(pkldict):
    #use analyze code
    section = 'soma0'
    V_section, t_section = np.array(pkldict[section]['Vm']), np.array(pkldict[section]['t'])
    t_stim = t_section[np.array(pkldict[section]['stimstate']>0)]

    #plt.plot(t_section*1e3,V_section)
    #plt.show()

    Vmax_index = argrelextrema(np.array(V_section), np.greater)[0]
    Vmax_t = t_section[Vmax_index]
    Vmax_V = V_section[Vmax_index]
    spiking_times = Vmax_t[Vmax_V>0]
    spiking_and_bounds = np.append(np.append(t_stim[0],spiking_times),t_stim[-1])
    ISI = np.diff(spiking_times)*1e3
    ISI_with_bounds = np.diff(spiking_and_bounds)*1e3
    nAPs = len(spiking_times)
    print('\n\nDEBUG:')
    for e in dir():
        print(e,sys.getsizeof(e))
    if nAPs > 0:
        return 1
    else:
        return 0

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


def main():
    #TODO: add these variables to parser
    # Parse command line arguments
    parser = AStimRealisticNeuronParser()
    args = parser.parse()
    args['method'] = [None]
    logger.setLevel(args['loglevel'])
    if args['mpi']:
        logger.warning('NEURON multiprocessing disabled')

    #START DEBUGGING VALUES (normally this should be given in the command line)
    fs_array = [0.75] #75%                                                                                         #default: 1 (100%)
    radius_array = args['radius'] #[64*1e-9] #[16*1e-9, 32*1e-9, 64*1e-9] # #16nm                                  #default: 3.2e-08 nm
    freq_array = args['freq'] # [100*1e3]#[100*1e3, 500*1e3, 3000*1e3] # #100kHz                                   #default: 500000. Hz
    args['section'] = ['myelin0',  'node0','soma0','unmyelin0','apical0','basal0', 'axon0'] # 7 type of sections   #default: None                                                                                         #default: None
    args['pltscheme'] = {'Vm': ['Vm'], 'Cm': ['Cm'], 'Qm': ['Qm'], 'iax' : ['iax']} #plotting variables            #default: None                                                                                     #default: 100000. Pa
    args['tstim'] = [0.1]                                                                                          #default: 0.0001 s
    args['toffset'] = [0.02]                                                                                       #default: 0.003 s
    #args['neuron'] = ['realneuron'] #this is actually not the way to change the neuron type
                                        #but this is irrelevant as 
    #args['nbursts'] = [2] #this argument needs to be changed for burst-mode                                       #default: 1
    DC_array = args['DC'] #[1.0] #[0.1,0.5, 1.] # #this argument needs to be changed for PW mode                   #default: 1.
    PRF_array = args['PRF'] # [100.] #[50., 100., 1000.] #

    #print(f'cmd arguments: \n{args}');quit()
    #END DEBUGGING VALUES


    # Run batch

    cell_nr = args['cell'][0]
    se = args['se'][0]
    stimulation = args['stimulation'][0] #"uniform" #"gaussian" #"gaussian3D"

    #with open('titrate','w') as ftit:
    #    ftit.write(f'Titration proces: \n\n\n')
    try:
        with open('titrate.pkl', 'rb') as fh:
            result_dict = pickle.load(fh)
    except:
        result_dict = {}
    #result_matr = np.zeros((len(fs_array), len(radius_array), len(freq_array), len(DC_array), len(PRF_array)))
    #result_dict = {'refs': {'fs': fs_array, 'radius': radius_array, 'freq': freq_array, 'DC': DC_array, 'PRF': PRF_array}, 'table': result_matr}


    logger.info('Starting fiber A-STIM simulation batch')
    #print(AStimParser.parseSimInputs(args)) #shows how the arguments are parsed in an array and not as keyword arguments (kwargs)
    output = []
    for fiber_class in args['type']:
        fiber = fiber_class(cell_nr, se) #BREAKPOINT
        for ita,a in enumerate(radius_array):
            fiber.a = a
            for itfs,fs in enumerate(fs_array):
                fiber.fs = fs
                for x0 in args['x0']:
                    for sigma in args['sigma']:
                        print(f"source parameters: x0 ={x0}, sigma = {sigma}") if 'gaussian' in stimulation else None

                        for itfreq,freq in enumerate(freq_array):
                            #for radius in [32*1e-9]:
                                #for fs in [0.75]:
                            for itDC,DC in enumerate(DC_array):
                                for itPRF,PRF in enumerate(PRF_array):
                                    args['fs'] = [fs] #75%                                                                         
                                    args['radius'] = [a] #16nm                                                                        
                                    args['freq'] = [freq] #100kHz                                                                           
                                    args['DC'] = [DC] #this argument needs to be changed for PW mode                           
                                    args['PRF'] = [PRF]

                                    #titration
                                    low_amp, high_amp = 0*1e3, 600*1e3
                                    amp = 10*1e3
                                    iter = 0
                                    #print(args)
                                    output_dict = simulate_realnrn(args,high_amp,stimulation,x0,sigma, fiber,iter)
                                    if not thresh_excited(output_dict):
                                        with open('titrate','a') as ftit:
                                            ftit.write(f'({freq}, {a}, {fs}, {DC}, {PRF}, {cell_nr}): unexcitable after {iter} iterations.\n\n')
                                        result_dict = add_dict_entry(result_dict, -10, [cell_nr, freq, a, fs, DC, PRF])
                                        continue
                                    while (low_amp == 0 or high_amp == 600*1e3):
                                        output_dict = simulate_realnrn(args,amp,stimulation,x0,sigma, fiber,iter)
                                        #with open('titrate','a') as ftit:
                                        #    ftit.write(f'Lower bound = {low_amp*1e-3}kPa, higher bound: {high_amp*1e-3}kPa, with a difference of: {(high_amp-low_amp)*1e-3}kPa, amplitude during the simulation was: {amp*1e-3}kPa\n')
                                        iter +=1
                                        if thresh_excited(output_dict):
                                            high_amp = amp
                                            amp = high_amp/2 
                                            if amp <= low_amp:
                                                break
                                            excited = 1
                                        else:
                                            low_amp = amp
                                            amp = 2*low_amp
                                            if amp >= high_amp:
                                                break
                                            excited = 0
                                        #with open('titrate','a') as ftit:
                                        #    ftit.write(f'Excited: {"True" if excited == 1 else "False"}. \nLower bound = {low_amp*1e-3}kPa, higher bound: {high_amp*1e-3}kPa, with a difference of: {(high_amp-low_amp)*1e-3}kPa, amp is now set to: {amp*1e-3}kPa\n\n')
                                        #print(f'Lower bound = {low_amp*1e-3}kPa, higher bound: {high_amp*1e-3}kPa, with a difference of: {(high_amp-low_amp)*1e-3}kPa, amp is now set to: {amp*1e-3}kPa')
                                    amp = (high_amp+low_amp)/2
                                    #with open('titrate','a') as ftit:
                                        #ftit.write(f'The calculated threshold amplitude after the binary search is: {amp*1e-3}kPa after {iter} iterations.\n\n\n')
                                    #print(f'The calculated threshold amplitude after the binary search is: {amp*1e-3}kPa after {iter} iterations.\n\n\n')

                                    epsilon = 100
                                    amp = (high_amp+low_amp)/2
                                    while (high_amp-low_amp) > epsilon:
                                        output_dict = simulate_realnrn(args,amp,stimulation,x0,sigma, fiber,iter)
                                        #with open('titrate','a') as ftit:
                                        #    ftit.write(f'Lower bound = {low_amp*1e-3}kPa, higher bound: {high_amp*1e-3}kPa, with a difference of: {(high_amp-low_amp)*1e-3}kPa, amplitude during the simulation was: {amp*1e-3}kPa\n')
                                        iter +=1
                                        if thresh_excited(output_dict):
                                            high_amp = amp
                                            excited = 1
                                        else:
                                            low_amp = amp
                                            excited = 0
                                        amp = (high_amp+low_amp)/2
                                        #with open('titrate','a') as ftit:
                                        #    ftit.write(f'Excited: {"True" if excited == 1 else "False"}. \nLower bound = {low_amp*1e-3}kPa, higher bound: {high_amp*1e-3}kPa, with a difference of: {(high_amp-low_amp)*1e-3}kPa, amp is now set to: {amp*1e-3}kPa\n\n')
                                        #print(f'Lower bound = {low_amp*1e-3}kPa, higher bound: {high_amp*1e-3}kPa, with a difference of: {(high_amp-low_amp)*1e-3}kPa, amp is now set to: {amp*1e-3}kPa')
                                    with open('titrate','a') as ftit:
                                        ftit.write(f'({freq}, {a}, {fs}, {DC}, {PRF}, {cell_nr}): {amp*1e-3}kPa after {iter} iterations.\n\n')
                                    result_dict = add_dict_entry(result_dict, amp, [cell_nr, freq, a, fs, DC, PRF])
                                    #result_dict['table'][itfs, ita, itfreq, itDC, itPRF] = amp
                                    #print(f'The calculated threshold amplitude after the binary search is: {amp*1e-3}kPa after {iter} iterations.\n\n')
    with open(r"C:\Users\jgazquez\RealSONIC\titrate.pkl", 'wb') as fh:
        pickle.dump(result_dict, fh)

if __name__ == '__main__': 
    main()
