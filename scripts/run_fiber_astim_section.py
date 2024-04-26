# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-09-27 14:28:52
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2023-03-29 18:39:34

''' Run simulations of an SENN SONIC fiber model with a specific point-neuron mechanism
    upon ultrasound stimulation at one node. '''

from PySONIC.core import Batch, NeuronalBilayerSonophore
from PySONIC.utils import logger
from PySONIC.parsers import AStimParser
from MorphoSONIC.core import SectionAcousticSource
from MorphoSONIC.parsers import SectionAStimFiberParser


def main():
    # Parse command line arguments
    parser = SectionAStimFiberParser()
    args = parser.parse()
    args['method'] = [None] #methods: full, hybrid or sonic
    logger.setLevel(args['loglevel']) #where does this argument come from?
    if args['mpi']: #multiprocessing
        logger.warning('NEURON multiprocessing disabled') #mp is enabled??

    # Run batch
    logger.info('Starting fiber A-STIM simulation batch')
    queue = [item[:2] for item in NeuronalBilayerSonophore.simQueue(
        *AStimParser.parseSimInputs(args), outputdir=args['outputdir'])] #create a (drive,protocol) queue
    if args['save']: #saves outputs when running simulation batches
        queue = [(item[0][:2], item[1]) for item in queue] #trim the queue so only (drive, protocol) and outputfolder is inside the queue
    output = []
    for fiber_class in args['type']: #'cortical_network', 'mrg', 'original_mrg', 'node_network', 'node_pop', 'ESTIM', 'radial_model', 'original_radial_model', 'senn', 'original_senn', 'smart_node_network', 'sweeney', 'original_sweeney', 'unmyelinated', 'original_unmyelinated'
        for fiberD in args['fiberD']: #fiber diameter
            for nnodes in args['nnodes']: #number of nodes
                fiber = fiber_class(fiberD, nnodes) #object instantiation
                for a in args['radius']: #sonophore radius
                    fiber.a = a
                    for fs in args['fs']: #sonophore coverage fraction
                        fiber.fs = fs
                        for sec_id in args['secid']: #ID of morphological section targeted by stimulus
                            if sec_id is None or sec_id == 'center':
                                sec_id = fiber.central_ID #only the central section is stimulated
                            if args['save']: #different 'formatting' of the items in the queue if the results need to be saved or not
                                simqueue = [(
                                    [SectionAcousticSource(
                                        sec_id, item[0][0].f, item[0][0].A), *item[0][1:]],
                                    item[1]
                                ) for item in queue]
                                func = fiber.simAndSave
                            else: 
                                simqueue = [
                                    [SectionAcousticSource(sec_id, item[0].f, item[0].A), *item[1:]]
                                    for item in queue]
                                func = fiber.simulate #this is the actual stimulation function
                            print(func)
                            batch = Batch(func, simqueue) #a batch is created with a simulate (and save) function and a simulation queue
                            output += batch(loglevel=args['loglevel']) #batch is added to the output
                            #quit()

    # Plot resulting profiles
    if args['plot'] is not None: # variables that need to be plotted -> 'all' plots all variables
        parser.parsePlot(args, output)


if __name__ == '__main__':
    main()
