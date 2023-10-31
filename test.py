""""just testing some things"""

import numpy as np
import re
import copy
import os
import sys
import argparse
# import tempFunctions as tf
# import PySONIC as ps

#     parser example
# def main():
#     # Create an ArgumentParser object
#     parser = argparse.ArgumentParser(description="Custom Argument Example")

#     # Add custom arguments for 't', 'f', and 'z' with default values
#     parser.add_argument("-t", type=int, default=1, help="Custom argument 't'")
#     parser.add_argument("-f", type=float, default=0.0, help="Custom argument 'f'")
#     parser.add_argument("-z", type=str, default="default_value", help="Custom argument 'z'")

#     # Parse the command-line arguments
#     args = parser.parse_args()

#     # Access the values of 't', 'f', and 'z' from the parsed arguments
#     t_value = args.t
#     f_value = args.f
#     z_value = args.z

#     # Now you can use these values in your script
#     print(f"t = {t_value}")
#     print(f"f = {f_value}")
#     print(f"z = {z_value}")

# if __name__ == "__main__":
#     main()

a = "/expanse/projects/nsg/home/nsguser/ngbw.expanse.nsgprod/workspace/NGBW-JOB-PY_CRI_SIMULATOR_EXPANSE-34C2320845C24688A2DE138C312D5A18/RealSIM/scripts"
b = a.split("RealSIM")[0]+"/RealSIM"
print(b)