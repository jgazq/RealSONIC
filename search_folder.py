import os

def search(string, folder):
    for root, dirs, files in os.walk(folder): #go through all files in the mechanics folders (all depths)
        for file in files:
            file = os.path.join(root,file)
            if file.endswith(".py") or file.endswith(".txt") or file.endswith(".hoc"):
                with open (file) as f:
                    for i,line in enumerate(f):
                        if string in line:
                            print(f'\" {string} \" found in \'{line}\' on line: {i+1}  -> in file: {file}\n')

# search('toPickle',os.getcwd())
# search('fromPickle',os.getcwd())
# search('insertVext',os.getcwd())
# search('getcwd()',os.getcwd())