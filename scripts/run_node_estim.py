# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2017-08-24 11:55:07
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2022-11-10 13:19:24

''' Run E-STIM simulations of a specific point-neuron. '''

from PySONIC.core import Batch, PointNeuron
from PySONIC.utils import logger
from PySONIC.parsers import EStimParser
from MorphoSONIC.models import Node


def main():
    # Parse command line arguments
    parser = EStimParser()
    args = parser.parse()
    logger.setLevel(args['loglevel'])
    if args['mpi']:
        logger.warning('NEURON multiprocessing disabled')
    sim_inputs = parser.parseSimInputs(args)
    simQueue_func = {5: 'simQueue', 6: 'simQueueBurst'}[len(sim_inputs)]

    # Run E-STIM batch
    logger.info("Starting node E-STIM simulation batch")
    queue = getattr(PointNeuron, simQueue_func)(
        *sim_inputs, outputdir=args['outputdir'], overwrite=args['overwrite'])
    output = []
    print(queue)
    for pneuron in args['neuron']:
        node = Node(pneuron)
        print(node)
        batch = Batch(node.simAndSave if args['save'] else node.simulate, queue)
        print(batch)
        output += batch(loglevel=args['loglevel'], mpi=args['mpi'])

    # Plot resulting profiles
    if args['plot'] is not None:
        parser.parsePlot(args, output)


if __name__ == '__main__':
    main()
