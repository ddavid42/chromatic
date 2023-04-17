#!/bin/python3
import argparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

parser = argparse.ArgumentParser(description='Read the cityplot')
parser.add_argument('lines', help='the number of lines')
parser.add_argument('columns', help='the number of columns')
parser.add_argument('file_data', help='the file name containing the data')
args = parser.parse_args()

fig = plt.figure()          #create a canvas, tell matplotlib it's 3d
ax = fig.add_subplot(111, projection='3d')

_x = np.arange(int(args.lines))
_y = np.arange(int(args.columns))
_xx, _yy = np.meshgrid(_x, _y)
x, y = _xx.ravel(), _yy.ravel()

bottom = np.zeros_like(x)
value = np.zeros_like(x)

i = 0
citymap = np.zeros((int(args.lines), int(args.columns)))
with open(args.file_data, "r") as data_file:
    for line in data_file:
        citymap[i] = [float(j) for j in line.split(",")][0:int(args.columns)]
        i = i+1
        if i >= int(args.lines):
            break

for i in range(int(args.lines)):
    for j in range(int(args.columns)):
        value[i * int(args.columns) + j] = citymap[i][j]

cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
max_height = np.max(value)   # get range of colorbars so we can normalize
min_height = np.min(value)
# scale each z to [0,1], and get their rgb values
rgba = [cmap((k-min_height)/max_height) for k in value] 

ax.bar3d(x, y, bottom, 1, 1, value, color=rgba, shade=True)
plt.title("chromatic analysis of the contributions of coefficients")
plt.xlabel("matrix width")
plt.ylabel("matrix length")
plt.show()
