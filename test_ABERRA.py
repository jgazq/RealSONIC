#from scipy.interpolate import griddata, interpn
#nrn_options = "-NSTACK 10000 -NFRAME 525"
#Aberra recommends: -NSTACK 100000 -NFRAME 20000
#os.environ["NEURON_MODULE_OPTIONS"] = nrn_options

import os
if "DISPLAY" in  os.environ:
    del os.environ['DISPLAY']

from neuron import h, gui
#h("NSTACK_size = 10000") 
"mosinit.command"

h.load_file('init.hoc')
#h.load_file('thresh4.hoc') 
#h.load_file('get_es2.hoc') #for TMS
#h.load_file('interp_coordinates.hoc') #for TMS (to calc Dx,Dy and Dz)


cell_nr = 7
h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
h.cell_chooser(cell_nr)

for sec in h.allsec():
    print(sec,sec.cm)