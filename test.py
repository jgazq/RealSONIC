""""just testing some things"""

import numpy as np
import re
import copy
# import tempFunctions as tf
# import PySONIC as ps


# with open('C:\\Users\\jgazquez\\RealSONIC\\mechanisms\\Ca.mod') as f:
#     text = f.read()
# print(text)
# search = re.search('PROCEDURE.*\{.*\}',text,flags=re.DOTALL)
# print(search)

d = {'a_': 5,'a_': 5, 'a_': 5, 'a_': 5 }
for key in d.keys():
    print(key)
    key = key.replace('_','')
    print(key)
print(d)