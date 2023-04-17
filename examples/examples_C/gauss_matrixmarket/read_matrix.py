#!/bin/python3
import argparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

parser = argparse.ArgumentParser(description='Read the matrix cityplot')
parser.add_argument('file_data', help='the file name containing the data')
parser.add_argument('lines', nargs='?', help='the number of lines')
parser.add_argument('columns', nargs='?', help='the number of columns')
args = parser.parse_args()

fig = plt.figure()          #create a canvas, tell matplotlib it's 3d
ax = fig.add_subplot(111, projection='3d')

citymap = None
lines = 0
columns = 0
with open(args.file_data, "r") as data_file:
    first_line = True
    i = 0
    for line in data_file:
        if line[0] == '%' and line[1] == '%':
            continue
        values = line.split()
        if first_line:
            lines = int(values[0])
            columns = int(values[1])
            count = int(values[2])
            first_line = False
            citymap = np.zeros((lines, columns))
        else:
            x = int(values[0])-1
            y = int(values[1])-1
            val = float(values[2])
            citymap[x][y] = val
            i = i+1
            if i >= count:
                break

nlines = args.lines
if not nlines:
    nlines = lines
else:
    nlines = int(nlines)

ncolumns = args.columns
if not ncolumns:
    ncolumns = columns
else:
    ncolumns = int(ncolumns)

_x = np.arange(nlines)
_y = np.arange(ncolumns)
_xx, _yy = np.meshgrid(_x, _y)
x, y = _xx.ravel(), _yy.ravel()

bottom = np.zeros_like(x)
value = np.zeros_like(x)

for i in range(lines):
    for j in range(columns):
        xpos = int((i + 0.5)*(float(nlines)/float(lines)))
        ypos = int((j + 0.5)*(float(ncolumns)/float(columns)))
        value[xpos * int(ncolumns) + ypos] += citymap[i][j]

cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
max_height = np.max(value)   # get range of colorbars so we can normalize
min_height = np.min(value)
# scale each z to [0,1], and get their rgb values
print('set color for citypot')
rgba = [cmap((k-min_height)/(max_height-min_height)) for k in value] 

print('set bars')
ax.bar3d(x, y, bottom, 1, 1, value, color=rgba, shade=True)
print('print figure')
plt.title("chromatic analysis of the contributions of coefficients")
plt.xlabel("matrix width")
plt.ylabel("matrix length")
plt.show()
