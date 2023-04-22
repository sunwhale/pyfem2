import os
import sys

sys.path.insert(0, os.getcwd())

import time

from pyfem.io.input import input_reader
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver import Solver

t1 = time.time()

props, globdat = input_reader()

# print(props)
# print(type(globdat))
#
solver = Solver(props, globdat)
output = OutputManager(props, globdat)

while globdat.active:
    solver.run(props, globdat)
    output.run(props, globdat)

t2 = time.time()

total = t2 - t1
print("Time elapsed = ", total, " [s].\n")

print("PyFem analysis terminated successfully.")
