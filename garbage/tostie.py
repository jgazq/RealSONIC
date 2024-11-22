import subprocess
import os
import sys

# Replace '1000' with the desired stack size
# Replace 'your_script.py' with the script you want to run in the new NEURON session

os.chdir('scripts')
subprocess.run('ls')
subprocess.run(["nrniv", "-python", "-NSTACK", "1000", "run_titrate_realistic_astim.py"])
