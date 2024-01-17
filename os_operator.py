import os
from docx import Document

def search_folder(string, folder,case_sens=1):
    """searches in all lines of all files in a folder for a given list of strings
    
    string: can be a list of different strings or a single string
    case_sens: the search is case sensitive by default
    """
    
    if type(string) != list:
        string = [string] 
    encounters = 0
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            if file.endswith(".py") or file.endswith(".txt") or file.endswith(".hoc"):
                if not 'os_operator' in file: #do not include the content of this file
                    with open (file) as f:
                        for i,line in enumerate(f):
                            if case_sens:
                                if all(word in line for word in string):
                                #if string in line:
                                    encounters += 1
                                    print(f'\nencounter {encounters}:\t\" {string} \" found in \'{line}\' on line: {i+1}  -> in file: {file}')
                            else: #case insensitive
                                if all(word.lower() in line.lower() for word in string):
                                #if string in line:
                                    encounters += 1
                                    print(f'\nencounter {encounters}:\t\" {string} \" found in \'{line}\' on line: {i+1}  -> in file: {file}')                               
    print('')

string_examples = ['toPickle', 'fromPickle', 'insertVext', 'getcwd()', 'criterion not met', '#to debug P_A = 0', 'MethodType', '#for RealDynNeuron', 'setMechValue', ['gating_from_PROCEDURES'], 'ABERRA', 'insert']
folder_examples = [os.getcwd(), r'C:\Users\jgazquez\PySONIC', r'C:\Users\jgazquez\MorphoSONIC']
search_folder(['listo'],r'C:\Users\jgazquez\MorphoSONIC')


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


def search_docx(string,folder):
    """searches in all lines of all word files in a folder for a given string"""

    encounters = 0
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            if file.endswith(".docx"):
                try:
                    f = Document(file)
                except:
                    print(file)
                    print('|__> this file is open in Word\n')
                for i,line in enumerate(f.paragraphs):
                    if string in line.text:
                        encounters += 1
                        print(f'\nencounter {encounters}:\t\" {string} \" found in \'{line.text}\' on line: {i+1}  -> in file: {file}')
    print('')

# search_docx('rest',r'C:\Users\jgazquez\OneDrive - UGent\PhD') #just a test
# search_docx('file is not open',r'C:\Users\jgazquez\OneDrive - UGent\PhD') #just a test