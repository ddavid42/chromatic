# Importing NumPy Library
import numpy as np
import matplotlib.patches as mpatches
import seaborn as sns               # For the automatic handling of gradient for the heatmap  
import sys
from random import *
from num2 import ChNbr
from math import *
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons

def P(x,y):
    return 9.*x*x*x*x-y*y*y*y+2.*y*y

N=2000
STEP = 5
# Making numpy array of N x N size and initializing 
# to zero for storing augmented matrix
a = np.zeros((N,N), dtype=np.dtype(ChNbr))

x = ChNbr(0.,follow=True)
y = ChNbr(0.,follow=True)

for i in range(0,N,STEP):
    for j in range(0,N,STEP):
        x = ChNbr(i, tab={1:1.0}, follow=False)
        y = ChNbr(j, tab={2:1.0}, follow=False)
        
        y.val = j
        a[i][j] = P(x,y) 

# Plot a 3D histogram of the results
print("Generate 3D plot")

# Creating dataset

X=np.zeros((ceil(N/STEP),ceil(N/STEP)))
Y=np.zeros((ceil(N/STEP),ceil(N/STEP)))
red=np.zeros((ceil(N/STEP),ceil(N/STEP)))
green=np.zeros((ceil(N/STEP),ceil(N/STEP)))
blue=np.zeros((ceil(N/STEP),ceil(N/STEP)))
val = np.zeros((ceil(N/STEP),ceil(N/STEP)))
cancellation = np.zeros((ceil(N/STEP),ceil(N/STEP)))

ii=0
for i in range(0,N,STEP):
    jj=0
    for j in range(0,N,STEP): 
        X[ii][jj]=i
        Y[ii][jj]=j
        v0 = a[i][j].tab[1]
        v1 = a[i][j].tab[2]
        v2 = a[i][j].tab[-1]
        red[ii][jj]=53*v0/(v0+v1+v2)
        green[ii][jj]=53*v1/(v0+v1+v2)
        blue[ii][jj]=53*v2/(v0+v1+v2)
        cancellation[ii][jj]= min(max(log2(fabs(v0+v1+v2)), 1) , 53)
        val[ii][jj]=a[i][j].val
        jj+=1
    ii+=1


# Creating figure
fig = plt.figure(figsize =(14, 9))
ax = plt.axes(projection ='3d')



# Adding labels
ax.set_xlabel('X-axis')
ax.set_xlim(0, N)
ax.set_ylabel('Y-axis')
ax.set_ylim(0, N)
ax.set_zlabel('Z-axis')
ax.set_zlim(np.min(red), np.max(red))
ax.set_title('Impact of x,y, constant and cancellation in the computation of F(x,y)=9.x**4 - y**4 + 2.y**2')

# Creating plot
s1 = ax.plot_surface(X, Y, red, visible=False, label = "x")
s1._facecolors2d=s1._facecolor3d
s1._edgecolors2d=s1._edgecolor3d
s2 = ax.plot_surface(X, Y, green, visible=False, label = "y")
s2._facecolors2d=s2._facecolor3d
s2._edgecolors2d=s2._edgecolor3d
s3 = ax.plot_surface(X, Y, blue, visible=False, label = "constant")
s3._facecolors2d=s3._facecolor3d
s3._edgecolors2d=s3._edgecolor3d
s4 = ax.plot_surface(X, Y, cancellation, visible=True, label = "cancellation")
s4._facecolors2d=s4._facecolor3d
s4._edgecolors2d=s4._edgecolor3d

surfs = [s1,s2,s3,s4]

plt.legend()

rax = fig.add_axes([0.05, 0.4, 0.1, 0.15])
labels = [str(surf.get_label()) for surf in surfs]
visibility = [surf.get_visible() for surf in surfs]
check = CheckButtons(rax, labels, visibility)

def func(label):
    index = labels.index(label)
    surfs[index].set_visible(not surfs[index].get_visible())
    plt.draw()

check.on_clicked(func)

# show plot
plt.show()
        