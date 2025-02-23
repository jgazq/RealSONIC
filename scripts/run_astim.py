# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2017-02-13 18:16:09
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2020-04-17 20:23:11

''' Run A-STIM simulations of a specific point-neuron. '''

from PySONIC.core import NeuronalBilayerSonophore, Batch
from PySONIC.utils import logger
from PySONIC.parsers import AStimParser, getPointNeuron


def main():
    # Parse command line arguments
    parser = AStimParser()
    args = parser.parse()
    print(args); #quit()
    args['neuron'] = [getPointNeuron('soma')]
    args['fs'] = [.75]
    args['radius'] = [64*1e-9] 
    args['freq'] = [100*1e3] 
    args['amp'] = [600*1e3] 
    args['tstim'] = [0.1]
    args['toffset'] = [0.01]
    args['plot'] = True   
    args['pltscheme'] = {'Vm': ['Vm'], 'Cm': ['Cm'], 'Qm': ['Qm']}
    logger.setLevel(args['loglevel'])
    sim_inputs = parser.parseSimInputs(args)
    simQueue_func = {9: 'simQueue', 10: 'simQueueBurst'}[len(sim_inputs)]

    # Run A-STIM batch
    logger.info("Starting A-STIM simulation batch")
    queue = getattr(NeuronalBilayerSonophore, simQueue_func)(
        *sim_inputs, outputdir=args['outputdir'], overwrite=args['overwrite'])
    output = []
    for a in args['radius']:
        for pneuron in args['neuron']:
            nbls = NeuronalBilayerSonophore(a, pneuron)
            batch = Batch(nbls.simAndSave if args['save'] else nbls.simulate, queue)
            output += batch(mpi=args['mpi'], loglevel=args['loglevel'])

    # Plot resulting profiles
    if args['plot'] is not None:
        parser.parsePlot(args, output)


if __name__ == '__main__':
    main()
