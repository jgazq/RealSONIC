""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt
import os

#this is added to work on macOS
if "DISPLAY" in  os.environ:
    del os.environ['DISPLAY']

from neuron import h
import tempFunctions as tf
import tempConstants as tc
