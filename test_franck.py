# Importing NumPy Library
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns               # For the automatic handling of gradient for the heatmap  
import sys
from random import *
from num2 import ChNbr

 
a = ChNbr(1, follow=True)

z = a + 100
r = z - 0.99*z

print(r)
