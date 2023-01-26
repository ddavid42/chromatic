# Importing NumPy Library
import sys
sys.path.append('../../..')
import numpy as np

from io import StringIO             # To import matrix in martrixmarket format
from scipy.io import mmread


import tracemalloc
import linecache
import time
import os
from collections import Counter

from random import *
from chnbr import ChNbr

def display_top(snapshot, key_type='lineno', limit=3):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

def Gauss( a, x , n):
    # Applying Gauss Elimination
    for i in range(n):
        print("GE(i): ",i)
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

tracemalloc.start()




# Reading number of unknowns
#n = int(input('Enter number of unknowns: '))
# EXEMPLE 1: of a 2x2 Ill conditionned system 
# n = 2
# b = np.array([[1.,1.,2.],[1.,1.001,2]])

#m = mmread('bcspwr01.mtx.gz')
m = mmread('bcsstk14.mtx.gz')
if (m.A.shape[0] == m.A.shape[1]):
    n = m.A.shape[1]
    print("Handling matrix of shape :", m.A.shape)
else:
    print("Not square matrix !")
    exit()

a = np.array(m.A)
b = np.array([i for i in range(n)])


start = time.perf_counter()
x = np.linalg.solve(a, b)
end = time.perf_counter()
print(f"\n#(Linalg) Time taken is {end - start}")

aa = np.zeros((n,n), dtype=np.dtype(ChNbr))
bb = np.zeros((n), dtype=np.dtype(ChNbr))
# Reading augmented matrix coefficients
for i in range(n):
    for j in range(n):
        aa[i][j] = ChNbr(a[i][j],follow=True)
        #aa[i][j] = a[i][j]

for i in range(n):
    #bb[i] = b[i]
    bb[i] = ChNbr(b[i],follow=True)

start = time.perf_counter()
xx = np.linalg.solve(aa, bb)
end = time.perf_counter()

print(f"\n#(Linalg 2) Time taken is {end - start}")

exit()

# Making numpy array of n x n+1 size and initializing 
# to zero for storing augmented matrix
#a = np.zeros((n,n+1), dtype=np.dtype(ChNbr))
a = np.zeros((n,n+1), dtype=np.dtype(float))

# Making numpy array of n size and initializing 
# to zero for storing solution vector
#x = np.zeros(n, dtype=np.dtype(ChNbr))
x = np.zeros(n, dtype=np.dtype(float))

# Reading augmented matrix coefficients
for i in range(n):
    for j in range(n):
        #a[i][j] = ChNbr(b[i][j],follow=True)
        a[i][j] = b[i][j]
    # last column B vector
    #a[i][n]=ChNbr(float(i+1),follow=True)
    a[i][n]=float(i+1)


#print("Initial System")
#print_np_mat_ChNbr(a)
 
Gauss( a, x, n )



end = time.perf_counter()

snapshot = tracemalloc.take_snapshot()
display_top(snapshot)

print("Time taken is ", (end - start))

