# Importing NumPy Library
from ast import Str
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import seaborn as sns               # For the automatic handling of gradient for the heatmap  
import sys
sys.path.append('../../..')

from random import *
from chnbr import ChNbr

MAXITER=200

# 1. Computation using float
# First constant to follow
golden_ratio = 1.+math.sqrt(5)/2
xi = 2.
xii = 0.

print("FLOAT COMPUTATION")
for i in range(0,MAXITER):
    xii = math.sqrt(xi+1)
    res =  (xi-golden_ratio)/golden_ratio
    print(i, xi, res)
    xi = xii

# 1. Computation using ChNbr
# First constant to follow
golden_ratio = ChNbr(1.+math.sqrt(5)/2, follow=True)
xi = ChNbr(2., follow=True)
xii = ChNbr(0., follow=True)

print("CHNBR COMPUTATION")
for i in range(0,MAXITER):
    xii = ChNbr.sqrt(xi+1)
    res =  (xi-golden_ratio)/golden_ratio
    print(i, xi.val, res.val)
    xi = xii
