import glob, os

example_drivers = glob.glob('*/*.py')

commands = []
try:
    os.remove('check_coverage.py')
except:
    pass

with open('check_coverage.py', 'w+') as fout:
    
    fout.write('import os\n')
    fout.write('import networkx as nx\n')
    fout.write('from pydfnworks import *\n\n\n')
    for name in example_drivers:
        dirname = name.split('/')
        # print(dirname)
        if dirname[0] == '4_user_rects':
            function = 'four_user_rests'
            fout.write(f'def {function}():\n')

        elif dirname[0] == '4_user_ell_uniform':
            function = 'four_user_ell_uniform'
            fout.write(f'def {function}():\n')
        else:
            fout.write(f'def {dirname[0]}():\n')
            commands.append(f"{dirname[0]}()")

        fout.write(f"\tprint('starting on {dirname[0]}')\n")

        # print(name)
        with open(name, 'r+') as fin:
            fout.write(f"\tsrc_path = '{dirname[0]}'\n")
            fout.write(f"\tos.chdir(src_path)\n")
            for line in fin.readlines():
                if not 'import' in line: # and not 'getcwd' in line:
                    fout.write(f"\t{line}")
        fout.write("\tos.chdir('/Users/jhyman/src/dfnWorks/examples/')")
        fout.write("\n\n")
        
    for c in commands:
        fout.write(f"{c}\n")

