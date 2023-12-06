import os

def search_folder(string, folder):
    """searches in all lines of all files in a folder for a given string"""
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            if file.endswith(".py") or file.endswith(".txt") or file.endswith(".hoc"):
                with open (file) as f:
                    for i,line in enumerate(f):
                        if string in line:
                            print(f'\" {string} \" found in \'{line}\' on line: {i+1}  -> in file: {file}\n')

# search_folder('toPickle',os.getcwd())
# search_folder('fromPickle',os.getcwd())
# search_folder('insertVext',os.getcwd())
# search_folder('getcwd()',os.getcwd())
# search_folder('criterion not met',os.getcwd())


def search_file(string, folder,filename_only=False):
    """searches in all files in a folder for a given string in the filename
    filename_only: string can only appear in the filename itself and not in the path"""
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            if not filename_only:
                file = os.path.join(root,file)
            if file.endswith(".py") or file.endswith(".txt") or file.endswith(".hoc"):
                if string in file:
                    print(os.path.join(root,file))

# search_file('neuron',os.getcwd(),1)


def remove_eff(folder):
    """remove the effective dulicates of mechanisms
    inverse operator of write_modl.py"""
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            if "eff" in file and file.endswith(".mod"):
                os.remove(file)

#remove_eff(os.path.join(os.getcwd(),'mechanisms'))


def remove_dll(folder):
    """remove the generated files to compile the mechanisms
    inverse operator of mknrndll"""
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            if file.endswith(".o") or file.endswith(".c") or file.endswith(".tmp") or file.endswith(".dll"):
                os.remove(file)
                
#remove_dll(os.path.join(os.getcwd(),'mechanisms\eff'))