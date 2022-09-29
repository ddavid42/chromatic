# Importing NumPy Library
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns               # For the automatic handling of gradient for the heatmap  
import sys
from random import *
from num2 import ChNbr

# Reading number of unknowns
#n = int(input('Enter number of unknowns: '))
# EXEMPLE 1: of a 2x2 Ill conditionned system 
# n = 2
# b = np.array([[1.,1.,2.],[1.,1.001,2]])

# EXEMPLE 2: of a 4x4 Ill conditionned system (Hilbert Matrix) 
n = 4
b = np.array( [[5., 7., 6., 5., 1.],
                [7., 10., 8., 7., 2.],
                [6., 8., 10., 9., 3.],
                [5., 7., 9., 10., 4.]])


# Making numpy array of n x n+1 size and initializing 
# to zero for storing augmented matrix
a = np.zeros((n,n+1), dtype=np.dtype(ChNbr))

# Making numpy array of n size and initializing 
# to zero for storing solution vector
x = np.zeros(n, dtype=np.dtype(ChNbr))

def print_np_mat_ChNbr(a):
    d = a.shape
    print("INDEX")
    for i in range(d[0]):
        for j in range(d[1]):
            print('%5d ' %(a[i][j].idx), end='')
        print("\n")

    print("VALUE")
    for i in range(d[0]):
        for j in range(d[1]):
            print(a[i][j].val, end='\t')
        print("\n")

# Reading augmented matrix coefficients
for i in range(n):
    for j in range(n+1):
        #a[i][j] = ChNbr( uniform(1,10) , follow=True)
        a[i][j] = ChNbr(b[i][j],follow=True)


print("Initial System")
print_np_mat_ChNbr(a)
    
# Applying Gauss Elimination
for i in range(n):
    if a[i][i] == 0.0:
        sys.exit('Divide by zero detected!')
        
    for j in range(i+1, n):
        ratio = a[j][i]/a[i][i]
        
        for k in range(n+1):
            a[j][k] = a[j][k] - ratio * a[i][k]


# Back Substitution
x[n-1] = a[n-1][n]/a[n-1][n-1]

for i in range(n-2,-1,-1):
    x[i] = a[i][n]
    
    for j in range(i+1,n):
        x[i] = x[i] - a[i][j]*x[j]
    
    x[i] = x[i]/a[i][i]

# Displaying solution
print('\nRequired solution is: ')
for i in range(n):
    print('X%d' %(i), x[i], "\n")

# Plot the Norm1 of each input values on the resulting system (ie. Detailed Condition Number)
heatmap1 = np.zeros((n,n+1))
for i in range(n):
    for k,v in x[i].tab.items():
        heatmap1[((k-1)//(n+1))][(k-1)%(n+1)]+=v

ax = sns.heatmap(heatmap1, annot=True, linewidth=0.5)
plt.show()