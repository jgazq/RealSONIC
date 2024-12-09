# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2017-06-02 17:50:10
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2023-04-06 10:16:57

''' Create lookup table for specific neuron. '''

import os
import itertools
import logging
import numpy as np
import sys
import datetime
sys.path.append(os.path.abspath(\
os.path.split(__file__)[0].split('RealSONIC')[0]+'\\RealSONIC'))
current_time = datetime.datetime.now()
now = datetime.datetime.strftime(current_time,'%Y_%m_%d_%H_%M_%S')

# from neuron import h, gui
# h.load_file('init.hoc')

from PySONIC.utils import logger, isIterable
from PySONIC.core import NeuronalBilayerSonophore, Batch, Lookup, AcousticDrive
from PySONIC.parsers import MechSimParser
from PySONIC.neurons import getDefaultPassiveNeuron
from PySONIC.constants import DQ_LOOKUP

# sys.path.append("c:\\users\\jgazquez\\RealSONIC")
# import tempFunctions as tf
# import tempConstants as tc
# import prev.functions as fs
# import prev.Interp3Dfield as tt


def computeAStimLookup(pneuron, aref, fref, Aref, Cm0ref, fsref, Qref, novertones=0,
                       test=False, mpi=False, loglevel=logging.INFO):
    ''' Run simulations of the mechanical system for a multiple combinations of
        imposed sonophore radius, US frequencies, acoustic amplitudes charge densities and
        (spatially-averaged) sonophore membrane coverage fractions, compute effective
        coefficients and store them in a dictionary of n-dimensional arrays.

        :param pneuron: point-neuron model
        :param aref: array of sonophore radii (m)
        :param fref: array of acoustic drive frequencies (Hz)
        :param Aref: array of acoustic drive amplitudes (Pa)
        :param Qref: array of membrane charge densities (C/m2)
        :param fsref: acoustic drive phase (rad)                                    #this should be sonophore coverage fraction?
        :param mpi: boolean statting wether or not to use multiprocessing
        :param loglevel: logging level
        :return: lookups dictionary
    '''
    descs = {
        'a': 'sonophore radii',
        'f': 'US frequencies',
        'A': 'US amplitudes',
        'fs': 'sonophore membrane coverage fractions',
        'overtones': 'charge Fourier overtones',
        'Cm0': 'membrane capacitance at rest'
    }

    # Populate reference vectors dictionary
    refs = {
        'a': aref,  # nm
        'f': fref,  # Hz
        'A': Aref,  # Pa
        'Q': Qref,  # C/m2
        'Cm0': Cm0ref # F/m2
    }

    err_span = 'cannot span {} for more than 1 {}'
    # If multiple sonophore coverage values, ensure that only 1 value of
    # sonophore radius and US frequency are provided
    if fsref.size > 1 or fsref.shape[0] != 1.: #checks if multiple fs are given, otherwise it needs to be 1
        for x in ['a', 'f']:
            assert refs[x].size == 1, err_span.format(descs['fs'], descs[x])
    # Add sonophore coverage vector to references
    refs['fs'] = fsref

    # If charge overtones are required, ensure that only 1 value of
    # sonophore radius, US frequency and coverage fraction are provided
    if novertones > 0:
        for x in ['a', 'f', 'fs']:
            assert refs[x].size == 1, err_span.format(descs['overtones'], descs[x])

    # If charge overtones are required, downsample charge and US amplitude input vectors
    if novertones > 0:
        nQmax = 50
        if len(refs['Q']) > nQmax:
            print("nQmax exceeded!\n\n")
            refs['Q'] = np.linspace(refs['Q'][0], refs['Q'][-1], nQmax) #same Q_m range but with fewer amount of samples: #50
        nAmax = 15
        if len(refs['A']) > nAmax:
            print("nAmax exceeded!\n\n")
            refs['A'] = np.insert(
                np.logspace(np.log10(refs['A'][1]), np.log10(refs['A'][-1]), num=nAmax - 1),
                0, 0.0) #Same P_A range but fewer amount of samples and conserving the 0 amplitude

    # If test case, reduce all vector dimensions to their instrinsic bounds
    if test:
        refs = {k: np.array([v.min(), v.max()]) if v.size > 1 else v for k, v in refs.items()} #in case of array: reduce to 2 values, otherwise keep the single value

    # Check validity of all reference vectors
    for key, values in refs.items():
        if not isIterable(values):
            raise TypeError(f'Invalid {descs[key]} (must be provided as list or numpy array)')
        if not all(isinstance(x, float) for x in values):
            raise TypeError(f'Invalid {descs[key]} (must all be float typed)')
        if len(values) == 0:
            raise ValueError(f'Empty {key} array')
        if key in ('a', 'f', 'Cm0') and min(values) <= 0:
            raise ValueError(f'Invalid {descs[key]} (must all be strictly positive)')
        if key in ('A', 'fs') and min(values) < 0:
            raise ValueError(f'Invalid {descs[key]} (must all be positive or null)')

    # Create simulation queue per sonophore radius
    drives = AcousticDrive.createQueue(refs['f'], refs['A']) #IMPORTANT LINE: creation of acoustic source(s), no phase given -> all possible combinations are created in the Batch class
    queue = []
    for drive in drives:
        for Qm in refs['Q']:
            #for Cm0 in refs['Cm0']: # we want to combine all the Cm0 values into 1 queue element so it can be combined into 1 calculation for Z_cycle
            queue.append([drive, refs['Cm0'], refs['fs'], Qm]) #queue contains all possible combinations of f, A (inside the drive), 'Cm0',fs and Qm

    # Add charge overtones to queue if required
    if novertones > 0:
        # Default references
        nAQ, nphiQ = 5, 5
        AQ_ref = np.linspace(0, 100e-5, nAQ)  # C/m2 #membrane charge amplitude
        phiQ_ref = np.linspace(0, 2 * np.pi, nphiQ, endpoint=False)  # rad #membrane charge phase
        # Downsample if test mode is on
        if test:
            AQ_ref = np.array([AQ_ref.min(), AQ_ref.max()]) #reduce to 2 values
            phiQ_ref = np.array([phiQ_ref.min(), phiQ_ref.max()]) #reduce to 2 values
        # Construct refs dict specific to Qm overtones
        Qovertones_refs = {}
        for i in range(novertones):
            Qovertones_refs[f'AQ{i + 1}'] = AQ_ref #assign the same array to every overtone
            Qovertones_refs[f'phiQ{i + 1}'] = phiQ_ref #assign the same array to every overtone
        # Create associated batch queue
        Qovertones = Batch.createQueue(*Qovertones_refs.values()) #creates all possible combinations of overtones values (each overtone has 5 possible A's and phi's) so 5^(2*novertones)
        Qovertones = [list(zip(x, x[1:]))[::2] for x in Qovertones] #this line converts every even set of variables in tuples of 2: [a,b,c,d] -> [(a,b),(c,d)]
        # Merge with main queue (moving Qm overtones into kwargs)
        print(len(queue),len(Qovertones))
        queue = list(itertools.product(queue, Qovertones)) #creates the "product" of two lists: ((x,y) for x in A for y in B)
        queue = [(x[0], {'Qm_overtones': x[1]}) for x in queue] #adding some encapsulation two the second element (overtones)
        # Update main refs dict, and reset 'fs' as last dictionary key
        refs.update(Qovertones_refs) #the overtone refs are added to the refs
        refs['Cm0'] = refs.pop('Cm0') #this puts Cm0 at the end (before fs), this is done to handle with the reshaping of tcomps
        refs['fs'] = refs.pop('fs') #this doesn't change anyting except the order of the dictionary, puts fs at the end (behind Qovertones_refs)

    # Get references dimensions
    dims = np.array([x.size for x in refs.values()]) #dimensions of ['a', 'f', 'A', 'Q', "Cm0", 'AQ1', 'phiQ1', ..., 'AQn', 'phiQn', 'fs']
    #tcompdims = np.array([refs[x].size for x in refs if (x != 'Cm0' and x != 'fs')])

    # Print queue (or reduced view of it)
    logger.info('batch queue:')
    Batch.printQueue(queue)

    # Run simulations and populate outputs
    logger.info('Starting simulation batch for %s neuron', pneuron.name)
    outputs = []
    for a in refs['a']:
        if pneuron.is_passive:
            xfunc = lambda *args, **kwargs: NeuronalBilayerSonophore(
                a, getDefaultPassiveNeuron()).computeEffVars(*args, **kwargs) #list with computation time and a list of dictionaries of effective variables
            batch = Batch(xfunc, queue)
        else:
            nbls = NeuronalBilayerSonophore(a, pneuron)
            batch = Batch(nbls.computeEffVars, queue)
        outputs += batch(mpi=mpi, loglevel=loglevel) #effective variables are calculated and then added to output together with the computation times
    # Split comp times and effvars from outputs
    effvars, tcomps = [list(x) for x in zip(*outputs)] #put the effvars dicts in a list, put the tcomp times in a list -> both effvars as tcomps is a list with length the number of possible combinations calculated
                                                       #length = #a x #f x #A x #Q x #AQ1 x #phiQ1 x ... x #AQn x #phiQn x #fs
    effvars = list(itertools.chain.from_iterable(effvars)) #effvars is a list that contains multiple lists of which the length is 1: the effvars dictionary -> simplify this to a list of dictionaries

    # Make sure outputs size matches inputs dimensions product
    nout = len(effvars) #output dimensions
    ncombs = dims.prod() #input dimensions
    if nout != ncombs:
        raise ValueError(
            f'Number of outputs ({nout}) does not match number of input combinations ({ncombs})')

    # Reshape effvars into nD arrays and add them to lookups dictionary
    logger.info(f'Reshaping {nout}-entries output into {tuple(dims)} lookup tables')
    varkeys = list(effvars[0].keys()) #list containing all the calculated output variables (this is the same for all input combinations so we take it from the first one)
    tables = {}
    for key in varkeys:
        effvar = [effvars[i][key] for i in range(nout)]
        tables[key] = np.array(effvar).reshape(dims)
    #tables is a dictionary in which each varkey(output variable) has the value of a multi-dimensional(dim = #inputs) list containing all the results for the different input variables

    # Reshape computation times, tile over extra fs dimension, and add it as a lookup table
    #DEBUG
    # print(refs.keys(),refs.values(),dims,tcompdims,dims[:-1])
    # tcomps1 = np.array(tcomps).reshape(dims[:-1]) #Cm dim is still here
    # print(f"tcomps1: {tcomps1}\ntcomps1.size: {tcomps1.shape}")
    # tcomps2 = np.array(tcomps).reshape(tcompdims) #Cm dim is also removed
    # print(f"tcomps2: {tcomps2}\ntcomp2s.size: {tcomps2.shape}")
    # tcompsrearranged = np.moveaxis(np.array([tcomps for i in range(dims[-1])]), 0, -1) #an extra dimsension is added so every value is now a list itself
    # print(f"tcomps rearranged: {tcompsrearranged}\ntcomp2s rearranged.size: {tcompsrearranged.shape}")
    # quit()
    #END DEBUG
    
    tcomps = np.array(tcomps).reshape(dims[:-2]) #the tcomps list is rearranged so every input argument can be indexed (# input variables = # dimensions list), also the last dimension (fs) is removed
                                                 #this assumes that Cm0 and fs are at the end, alternatively tcompdims can be used but order is then unknown for further use
    #tcomps = np.moveaxis(np.array([tcomps for i in range(dims[-1])]), 0, -1) #an extra dimsension is added so every value is now a list itself (this is done for fs dimension)
    tcomps = np.moveaxis([[tcomps for i in range(dims[-1])] for i in range(dims[-2])],[0,1],[-2,-1]) #Cm is last added, needs to be at the penultimate place and fs at the end
    #tcomps = np.moveaxis(np.array([tcomps for i in range(dims[-1])]), 0, -1)
    #|-> this is actually not correct because tcomp takes total time over all Cm0 and fs values and copies this value for every Cm0 or fs value
    #
    tables['tcomp'] = tcomps

    # Construct and return lookup object
    return Lookup(refs, tables)


def main():

    #print(f'{"^^"*284}\n\n{"^^"*284}\n\n')

    # cell_nr = 7
    # h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
    # h.cell_chooser(cell_nr)

    # cell_folder = "L23_PC_cADpyr229_2" #"L4_LBC_cACint209_2" #input("Give folder name of BBP cell:")
    # mech_folder = "cells/"+cell_folder+"/mechanisms/"
    # mod_files, mod_names = tf.read_mod(mech_folder)
    # l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names)

    parser = MechSimParser(outputdir='.') #why put the outputdir = .???
    parser.addNeuron()
    parser.addTest()
    parser.defaults['neuron'] = 'realneuron' #'RS' 
    parser.defaults['radius'] = np.array([16.0, 32.0, 64.0])  # nm #3 different sonophore radii
    parser.defaults['freq'] = np.array([20., 100., 500., 1e3, 2e3, 3e3, 4e3])  # kHz #7 frequencies
    parser.defaults['amp'] = np.insert(
        np.logspace(np.log10(0.1), np.log10(600), num=50), 0, 0.0)  # kPa #50 different pressure amplitudes, including 0 so none
    parser.defaults['charge'] = np.nan
    parser.defaults['Qstart'], parser.defaults['Qend'] = 0, -1
    parser.defaults['Cm0'] = np.array([1., 2.]) #np.array([0.02, 1., 2.]) #uF/cm2 #given in command line
    parser.add_argument('--novertones', type=int, default=0, help='Number of Fourier overtones') # numer of Qm overtones
    args = parser.parse()
    logger.setLevel(args['loglevel'])
    #print(args)

    for pneuron in args['neuron']: 

        #Give pneuron worst case Cm0 -> not anymore cuz calculations go wrong otherwise
        pneuron.Cm0 = max(args['Cm0']) #max(args['Cm0']) #this will give an error if 2 Cm0-values are given as input, only 1 can be given

        # Determine charge vector
        charges = args['charge']
        if charges.size == 1 and np.isnan(charges[0]):
            Qmin, Qmax = pneuron.Qbounds #np.array([np.round(self.Vm0 - 35.0), 50.0]) * self.Cm0 * 1e-3  # C/m2
            charges = np.arange(Qmin, Qmax + DQ_LOOKUP, DQ_LOOKUP)  # C/m2 #DQ_LOOKUP is step value
            #charges = np.linspace(Qmin, Qmax, 2) #for debugging
            Qstart, Qend = args['Qstart'][0], args['Qend'][0]
            if Qstart > 0:
                if Qend > 0:
                    charges = charges[Qstart:Qend]
                else:
                    charges = charges[Qstart:]
            elif Qend > 0:
                charges = charges = charges[:Qend]
        # Number of Fourier overtones
        novertones = args['novertones']

        # Determine output filename
        input_args = {'a': args['radius'], 'f': args['freq'], 'A': args['amp'], 'fs': args['fs'], 'Cm0': args['Cm0']}
        input_args = {**input_args, 'Qstart': args['Qstart'], 'Qend': args['Qend']} #addition of Qm0-range in the LUT name
        fname_args = {k: v[0] if v.size == 1 else None for k, v in input_args.items()}
        fname_args['novertones'] = novertones
        lookup_fpath = NeuronalBilayerSonophore(32e-9, pneuron).getLookupFilePath(**fname_args)
        fcode, fext = os.path.splitext(lookup_fpath)
        lookup_fpath = f'{fcode}_{now}{fext}' #add time to make each lookup calculation unique

        # Combine inputs into single list
        inputs = [args[x] for x in ['radius', 'freq', 'amp', 'Cm0', 'fs']] + [charges]

        # Adapt inputs and output filename if test case
        if args['test']:
            fcode, fext = os.path.splitext(lookup_fpath)
            lookup_fpath = f'{fcode}_test{fext}'

        # Check if lookup file already exists -> this becomes useless as timestamp of computations is added to the LUT .pkl file name
        if os.path.isfile(lookup_fpath):
            logger.warning(
                f'"{lookup_fpath}" file already exists and will be overwritten. Continue? (y/n)')
            user_str = input()
            if user_str not in ['y', 'Y']:
                logger.error('%s Lookup creation canceled', pneuron.name)
                return

        # Compute lookup
        lkp = computeAStimLookup(pneuron, *inputs, novertones=novertones,
                                 test=args['test'], mpi=args['mpi'], loglevel=args['loglevel']) #function defined above that actually computes the lookup file
        logger.info(f'Generated lookup: {lkp}')

        # Save lookup in PKL file
        logger.info('Saving %s neuron lookup in file: "%s"', pneuron.name, lookup_fpath)
        # print(f"V = {lkp.tables['V']}")
        # print(f"C = {lkp.refs['Q'] / lkp.tables['V']*1e5}")
        lkp.toPickle(lookup_fpath)


if __name__ == '__main__':
    main()
