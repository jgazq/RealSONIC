import os

def search(string, folder):
    for root, dirs, files in os.walk(folder):
        for file in files:
            with open (root+file,errors="ignore") as f:
                for line in f:
                    if string in line:
                        print(file,line)

search('get_axon_biophys()','C:\\Users\\jgazquez\\RealSONIC\\')