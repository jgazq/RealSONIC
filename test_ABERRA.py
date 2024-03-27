#from scipy.interpolate import griddata, interpn
#nrn_options = "-NSTACK 10000 -NFRAME 525"
#Aberra recommends: -NSTACK 100000 -NFRAME 20000
#os.environ["NEURON_MODULE_OPTIONS"] = nrn_options

from neuron import h, gui
#h("NSTACK_size = 10000") 
"mosinit.command"

h.load_file('init.hoc')
#h.load_file('thresh4.hoc') 
#h.load_file('get_es2.hoc') #for TMS
#h.load_file('interp_coordinates.hoc') #for TMS (to calc Dx,Dy and Dz)


cell_nr = 11
h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
h.cell_chooser(cell_nr)

# for sec in h.allsec():
#     print(sec,sec.cm)
reversals = {"ek": [], "ehcn": [], "ena": [], "eca": [], "e_pas": []}

for e in h.allsec():
    for rev, rev_vals in reversals.items():
        #print(e,rev)
        try:
            rev_val = eval(f"h.{e}.{rev}")
            if rev_val not in rev_vals:
                if len(rev_vals) != 0:
                    print("multiple rev_vals for the same reversal potential?")
                rev_vals.append(rev_val)
        except:
            pass

print(reversals)