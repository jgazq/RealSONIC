# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2017-06-13 09:40:02
# @Last Modified by:   Theo Lemaire
# @Last Modified time: 2023-05-18 09:06:45

import os
from setuptools import setup

readme_file = 'README.md'
requirements_file = 'requirements.txt'

def description():
    with open(readme_file, encoding="utf8") as f:
        started = False
        lines = []
        for line in f:
            if not started:
                if line.startswith('# Description'):
                    started = True
            else:
                if line.startswith('#'):
                    break
                else:
                    lines.append(line)
    desc = ''.join(lines).strip('\n')
    more = 'More information on the [https://github.com/tjjlemaire/PySONIC](Github repository).'
    return f'{desc}\n\n{more}'


def getFiles(path):
    return [f'{path}/{x}' for x in os.listdir(path)]


def get_requirements():
    with open(requirements_file, 'r', encoding="utf8") as f:
        reqs = f.read().splitlines() 
    return reqs


setup(
    name='PySONIC',
    version='1.0',
    description='Python implementation of the NICE & SONIC models, allowing for rapid simulations of ultrasound-evoked membrane dynamics in various point-neuron models.',
    long_description=description(),
    url='https://github.com/tjjlemaire/PySONIC',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    keywords=('SONIC NICE acoustic ultrasound ultrasonic neuromodulation neurostimulation excitation\
             computational model intramembrane cavitation'),
    author='Theo Lemaire',
    author_email='theo.lemaire1@gmail.com',
    license='MIT',
    packages=['PySONIC'],
    scripts=getFiles('scripts') + getFiles('tests'),
    install_requires=get_requirements(),
    zip_safe=False
)
