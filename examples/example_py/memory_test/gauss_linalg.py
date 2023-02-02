
# Importing NumPy Library
import numpy as np
import sys
import time

# Reading number of unknowns
n = 400

start = time.perf_counter()
a = np.array([[1/(i+j) for i in range(1,n+1)] for j in range(n)])
b = np.array([i for i in range(n)])
x = np.linalg.solve(a, b)
end = time.perf_counter()
print(f"\n#(Linalg) Time taken is {end - start}")


start = time.perf_counter()


# Making numpy array of n x n+1 size and initializing 
# to zero for storing augmented matrix
a = np.zeros((n,n+1))

# Making numpy array of n size and initializing 
# to zero for storing solution vector
x = np.zeros(n)

# Reading augmented matrix coefficients
a = np.array([[1/(i+j) for i in range(1,n+2)] for j in range(n)])

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
#print('\nRequired solution is: ')
#for i in range(n):
#    print('X%d = %0.2f' %(i,x[i]), end = '\t')



end = time.perf_counter()
print(f"\n#Time taken is {end - start}")

