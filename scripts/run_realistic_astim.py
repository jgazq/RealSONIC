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


def main():
    #TODO: add these variables to parser
    cell_nr = 7
    se = 0
    stimulation = "uniform"#"gaussian3D" #"gaussian"
    # Parse command line arguments
    parser = AStimRealisticNeuronParser()
    args = parser.parse()
    args['method'] = [None]
    logger.setLevel(args['loglevel'])
    if args['mpi']:
        logger.warning('NEURON multiprocessing disabled')

    #START DEBUGGING VALUES (normally this should be given in the command line)
    args['fs'] = [0.75]
    #END DEBUGGING VALUES

    # Run batch
    logger.info('Starting fiber A-STIM simulation batch')
    queue = [item[:2] for item in NeuronalBilayerSonophore.simQueue(
        *AStimParser.parseSimInputs(args), outputdir=args['outputdir'])]
    if args['save']:
        queue = [(item[0][:2], item[1]) for item in queue]
    output = []
    for fiber_class in args['type']:
        fiber = fiber_class(cell_nr, se)
        for a in args['radius']:
            fiber.a = a
            for fs in args['fs']:
                fiber.fs = fs
                for x0 in args['x0']:
                    for sigma in args['sigma']:
                        print(f"x0 ={x0}, sigma = {sigma}")
                        if args['save']:
                            if stimulation == "gaussian":
                                simqueue = [(
                                    [GaussianAcousticSource(
                                        x0, sigma, item[0][0].f, item[0][0].A), *item[0][1:]],
                                    item[1]
                                ) for item in queue]
                            elif stimulation == "uniform":
                                simqueue = [(
                                    [UniformAcousticSource(
                                        item[0][0].f, item[0][0].A), *item[0][1:]],
                                    item[1]
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
                        output += batch(loglevel=args['loglevel'])
    print(args['section'])
    refsec = fiber.sections[args['section'][0][:-1]][args['section'][0]] #fiber.refsection
    amplitude = gaussian(refsec.x_xtra,x0,sigma,queue[0].A) if stimulation == 'gaussian' else gaussian3D(refsec.x_xtra, refsec.y_xtra, refsec.y_xtra,x0,x0,x0,sigma,sigma,sigma,queue[0].A) if stimulation == 'gaussian3D' else queue[0][0].A
    print(f'{"-"*50}\nrefsection:\nlocation:\t({refsec.x_xtra}, {refsec.y_xtra}, {refsec.z_xtra})\namplitude = {amplitude}')
    args['plot'] = 'Vm' #for debugging the plot section
    # Plot resulting profiles
    if args['plot'] is not None:
        parser.parsePlot(args, output)


if __name__ == '__main__': 
    main()
