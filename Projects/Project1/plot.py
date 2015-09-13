# FILENAME: plot.py
import numpy as np
import sys
from Numeric import *
import Gnuplot
g = Gnuplot.Gnuplot(persist=1)
try:
infilename = sys.argv[1]
except:
print "Usage of this script", sys.argv[0], "infile", sys.argv[1]; sys.exit(1)
# Read file with data
ifile = open(infilename, 'r')
# Fill in x and y
x = [] ; y = []
for line in ifile:
pair = line.split()
x = float(pair[0]); y = float(pair[1])
ifile.close()
# convert to a form that the gnuplot interface can deal with
d = Gnuplot.Data(x, y, title='data from output file', with='lp')
g.xlabel('log10(h)') # make x label
g.ylabel('log10(|Exact-Computed|)/|Exact|')
g.plot(d)
# plot the data
g.hardcopy(filename="relerror.ps",terminal="postscript", enhanced=1, color=1)
