#!/bin/python3
import argparse
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns               # For the automatic handling of gradient for the heatmap  

parser = argparse.ArgumentParser(description='Read the headmap')
parser.add_argument('lines', help='the number of lines')
parser.add_argument('columns', help='the number of columns')
parser.add_argument('file_data', help='the file name containing the head data')
args = parser.parse_args()

i = 0
heatmap = np.zeros((int(args.lines), int(args.columns)))
with open(args.file_data, "r") as data_file:
    for line in data_file:
        heatmap[i] = [float(i) for i in line.split(",")]
        i = i+1
ax = sns.heatmap(heatmap, annot=True, linewidth=0.5)
plt.show()
