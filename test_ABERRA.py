#from scipy.interpolate import griddata, interpn
# nrn_options = "-NSTACK 100000 -NFRAME 525"
# # Aberra recommends: -NSTACK 100000 -NFRAME 20000
# import os
# os.environ["NEURON_MODULE_OPTIONS"] = nrn_options

from neuron import h, gui
import time
# h("NSTACK_size = 100000") 
"mosinit.command"

h.load_file('init.hoc')
#h.load_file('thresh4.hoc') 
#h.load_file('get_es2.hoc') #for TMS
#h.load_file('interp_coordinates.hoc') #for TMS (to calc Dx,Dy and Dz)


cell_nr = 2
h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
h.cell_chooser(cell_nr)
#time.sleep(200); quit() #sleep so  the NEURON interface can be used

for sec in h.allsec():
    print(sec,sec.nseg)
    # print(sec,sec.cm)
    #print(sec.psection()['species'])
    for seg in sec.allseg():
        print(seg)
    quit()

a = h.finitialize(-75)
print(a)
print('done')