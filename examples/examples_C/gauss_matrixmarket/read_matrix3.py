#!/bin/python3
import argparse
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

parser = argparse.ArgumentParser(description='Read the matrix cityplot')
parser.add_argument('file_data', help='the file name containing the data')
args = parser.parse_args()

fig = plt.figure()          #create a canvas, tell matplotlib it's 3d
ax = fig.add_subplot(111, projection='3d')

lines = 0
columns = 0
data = {}

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
        else:
            x = int(values[0])
            y = int(values[1])
            width = int(values[2])
            length = int(values[3])
            val = abs(float(values[4]))
            data[x * columns + y] = (width, length, val)
            data[y * columns + x] = (width, length, val)
            i = i+1
            if i >= count:
                break

_x = np.arange(lines)
_y = np.arange(columns)
_xx, _yy = np.meshgrid(_x, _y)
# x, y = _xx.ravel(), _yy.ravel()

x = np.zeros(len(data))
y = np.zeros(len(data))
width = np.zeros(len(data))
length = np.zeros(len(data))
bottom = np.zeros(len(data))
value = np.zeros(len(data))

i = 0
for k, v in data.items():
    x[i] = k / columns
    y[i] = k % columns
    width[i] = v[0]
    length[i] = v[1]
    value[i] = v[2]
    i = i+1

cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
max_height = np.max(value)   # get range of colorbars so we can normalize
min_height = np.min(value)
# scale each z to [0,1], and get their rgb data
print('set color for citypot')
rgba = [cmap((k-min_height)/(max_height-min_height)) for k in value] 

print('set bars')
ax.bar3d(x, y, bottom, width, length, value, color=rgba, shade=True)
print('print figure')
plt.title("chromatic analysis of the contributions of coefficients")
plt.xlabel("matrix width")
plt.ylabel("matrix length")
plt.show()
