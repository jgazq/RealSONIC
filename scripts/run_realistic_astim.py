# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-09-27 14:28:52
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2023-03-29 18:39:34

''' Run simulations of a multi-compartmental realistic Aberra/BBP neuron model with a specific point-neuron mechanism
    upon ultrasound stimulation at one node. '''

from PySONIC.core import Batch, NeuronalBilayerSonophore
from PySONIC.utils import logger
from PySONIC.parsers import AStimParser
from MorphoSONIC.core import GaussianAcousticSource
from MorphoSONIC.parsers import AStimRealisticNeuronParser


def main():
    #TODO: add these variables to parser
    cell_nr = 7
    se = 0
    # Parse command line arguments
    parser = AStimRealisticNeuronParser()
    args = parser.parse()
    args['method'] = [None]
    logger.setLevel(args['loglevel'])
    if args['mpi']:
        logger.warning('NEURON multiprocessing disabled')

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
                        if args['save']:
                            simqueue = [(
                                [GaussianAcousticSource(
                                    x0, sigma, item[0][0].f, item[0][0].A), *item[0][1:]],
                                item[1]
                            ) for item in queue]
                            func = fiber.simAndSave
                        else:
                            simqueue = [
                                [GaussianAcousticSource(x0, sigma, item[0].f, item[0].A), *item[1:]]
                                for item in queue]
                            func = fiber.simulate
                        batch = Batch(func, simqueue)
                        output += batch(loglevel=args['loglevel'])

    # Plot resulting profiles
    if args['plot'] is not None:
        parser.parsePlot(args, output)


if __name__ == '__main__': 
    main()
