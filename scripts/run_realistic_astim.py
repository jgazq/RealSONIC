# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-09-27 14:28:52
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2023-03-29 18:39:34

''' Run simulations of a multi-compartmental realistic Aberra/BBP neuron model with a specific point-neuron mechanism
    upon ultrasound stimulation at one node. '''

from PySONIC.core import Batch, NeuronalBilayerSonophore
from PySONIC.utils import logger, gaussian, gaussian3D
from PySONIC.parsers import AStimParser
from MorphoSONIC.core import GaussianAcousticSource, UniformAcousticSource, Gaussian3DAcousticSource
from MorphoSONIC.parsers import AStimRealisticNeuronParser
import pickle


def main():
    #TODO: add these variables to parser
    cell_nr = 7
    se = 0
    stimulation = "uniform" #"gaussian" #"gaussian3D"
    # Parse command line arguments
    parser = AStimRealisticNeuronParser()
    args = parser.parse()
    args['method'] = [None]
    logger.setLevel(args['loglevel'])
    if args['mpi']:
        logger.warning('NEURON multiprocessing disabled')

    #START DEBUGGING VALUES (normally this should be given in the command line)
    args['fs'] = [0.75] #75%                                                                                       #default: 1 (100%)
    args['radius'] = [64*1e-9] #16nm                                                                               #default: 3.2e-08 nm
    args['freq'] = [100*1e3] #100kHz                                                                               #default: 500000. Hz
    args['section'] = ['myelin0',  'node0','soma0','unmyelin0','apical0','basal0', 'axon0'] # 7 type of sections   #default: None
    args['plot'] = False                                                                                            #default: None
    args['pltscheme'] = {'Vm': ['Vm'], 'Cm': ['Cm'], 'Qm': ['Qm'], 'iax' : ['iax']} #plotting variables            #default: None
    args['amp'] = [100*1e3]                                                                                        #default: 100000. Pa
    args['tstim'] = [0.1]                                                                                          #default: 0.0001 s
    args['toffset'] = [0.01]                                                                                       #default: 0.003 s
    #args['neuron'] = ['realneuron'] #this is actually not the way to change the neuron type
                                     #but this is irrelevant as 
    #args['nbursts'] = [2] #this argument needs to be changed for burst-mode                                        #default: 1
    #args['DC'] = [.5] #this argument needs to be changed for PW mode                                               #default: 1.
    #args['PRF'] = [200.]

    #print(f'cmd arguments: \n{args}');quit()
    #END DEBUGGING VALUES

    # Run batch
    logger.info('Starting fiber A-STIM simulation batch')
    #print(AStimParser.parseSimInputs(args)) #shows how the arguments are parsed in an array and not as keyword arguments (kwargs)
    queue = [item[:2] for item in NeuronalBilayerSonophore.simQueue(
        *AStimParser.parseSimInputs(args), outputdir=args['outputdir'])]
    if args['save']:
        queue = [(item[0][:2], item[1]) for item in queue]
    output = []
    for fiber_class in args['type']:
        fiber = fiber_class(cell_nr, se) #BREAKPOINT
        for a in args['radius']:
            fiber.a = a
            for fs in args['fs']:
                fiber.fs = fs
                for x0 in args['x0']:
                    for sigma in args['sigma']:
                        print(f"source parameters: x0 ={x0}, sigma = {sigma}") if 'gaussian' in stimulation else None
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
    print(args['section'])
    refsec = fiber.sections[args['section'][0][:-1]][args['section'][0]] #fiber.refsection
    amplitude = gaussian(refsec.x_xtra,x0,sigma,queue[0][0].A) if stimulation == 'gaussian' else gaussian3D(refsec.x_xtra, refsec.y_xtra, refsec.y_xtra,x0,x0,x0,sigma,sigma,sigma,queue[0][0].A) if stimulation == 'gaussian3D' else queue[0][0].A
    print(f'{"-"*50}\nrefsection:\nlocation:\t({refsec.x_xtra}, {refsec.y_xtra}, {refsec.z_xtra}) um\namplitude = {amplitude*1e-3}kPa')
    args['ref_loc'] = (refsec.x_xtra, refsec.y_xtra, refsec.z_xtra)

    tosave = output[0][0].data
    outdir = r"C:\Users\jgazquez\RealSONIC\pkldump\dump_" + f"{args['fs'][0]*100}%_{args['radius'][0]*1e9}nm_{args['freq'][0]*1e-3}kHz" + \
        f"_{args['amp'][0]*1e-3}kPa_{args['tstim'][0]*1e3}ms_{args['toffset'][0]*1e3}ms_{args['PRF'][0]}Hz_{args['DC'][0]}DC" + ".pkl"
    with open(outdir, 'wb') as fh:
        pickle.dump(tosave, fh)

    # Plot resulting profiles
    if args['plot'] is not None:
        args['fiber'] = fiber
        parser.parsePlot(args, output) #BREAKPOINT


if __name__ == '__main__': 
    main()
