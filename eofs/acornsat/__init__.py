import os
import sys
cwd = os.getcwd()
os.chdir('..')
prev = os.getcwd()
sys.path.append(prev)
os.chdir(cwd)
