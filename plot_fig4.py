import os
import matplotlib.pyplot as plt
import matplotlib
from PySONIC.utils import si_format

from PySONIC.plt import plotEffectiveVariables
from PySONIC.utils import logger
from PySONIC.neurons import getPointNeuron

#added by Joaquin
project_root = os.path.dirname(os.path.realpath('utils.py'))
dataroot = project_root
#added by Joaquin

# Matplotlib parameters
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'arial'


def subdirectory(name):
    ''' Return (and create if needed) a data sub-directory. '''
    if dataroot is None:
        raise ValueError('You must specify the data root directory')
    if not os.path.isdir(dataroot):
        raise ValueError('You must create the data root directory prior to running the code')
    subdir = os.path.join(dataroot, name)
    if not os.path.isdir(subdir):
        os.mkdir(subdir)
    return subdir


def codes(a, pneuron, Fdrive, PRF, tstim):
    ''' Return codes dictionary for a combination of parameters. '''
    return [
        pneuron.name,
        f'{si_format(a, space="")}m',
        f'{si_format(Fdrive, space="")}Hz',
        f'PRF{si_format(PRF, space="")}Hz',
        f'{si_format(tstim, space="")}s'
    ]


def saveFigsAsPDF(figs, figindex, pneur_name): #changed by Joa
    ''' Save figures as PDFs with a specific figure index. '''
    figdir = subdirectory('figs')
    figbase = f'fig{figindex:02}'
    for fname, fig in figs.items():
        fig.savefig(os.path.join(figdir, f'{figbase}{fname}_{pneur_name}.pdf'), transparent=True) #changed by Joa


def subthr(x):
    ''' Sub-threshold acoustic pressure (Pa). '''
    return x - 5e3


def suprathr(x):
    ''' Supra-threshold acoustic pressure (Pa). '''
    return x + 20e3


figindex = 4
# fs = 12
# lw = 2
# ps = 15
figs = {}


pneur_name = 'realneuron'
pneuron = getPointNeuron(pneur_name)
a = 32e-9  # m
Fdrive = 500e3  # Hz
Adrive = 50e3  # Pa

fig = plotEffectiveVariables(pneuron, a=a, f=Fdrive, cmap='Oranges', zscale='log')
figs['a'] = fig

fig = plotEffectiveVariables(pneuron, a=a, A=Adrive, cmap='Greens', zscale='log')
figs['b'] = fig

fig = plotEffectiveVariables(pneuron, f=Fdrive, A=Adrive, cmap='Blues', zscale='log')
figs['c'] = fig

plt.show()
saveFigsAsPDF(figs, figindex,pneur_name)