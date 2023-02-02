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

# Number of iteration
N = 27

# 1. AS FLOAT
# Archimedes methods to compute PI
t1 =1/math.sqrt(3)

print("AS FLOAT")
for i in range(1, N):
    print(i,t1)
    t2 = (math.sqrt(t1*t1+1.0) - 1.0)/t1
    t1 = t2



# 1. AS CHNBR
# lst : List of chromatic numbers
lst=[]

print("AS CHNBR")
# Archimedes methods to compute PI
t1 =ChNbr(1/math.sqrt(3), follow=True)

for i in range(1, N):
    lst.append(t1)
    t2 = (ChNbr.sqrt(t1*t1+1.0) - 1.0)/t1
    t1 = t2

# Generate animated plot describing the evolution of t1 weight in the computed results
# at each iteration 


# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
ims = []

key_max1 = 1
key_max2 = -1

for it in range(0,N-1):
    # list out keys and values separately
    # key_list = list(lst[it].tab.keys())
    # val_list = list(lst[it].tab.values())

    # sorted_val_list = sorted(val_list)

    # key_max1 = key_list[val_list.index(sorted_val_list[-1])]
    # if len(val_list)>1:
    #     key_max2 = key_list[val_list.index(sorted_val_list[-2])]
    # else:
    #     key_max2 = 0

    # if len(val_list)>2:
    #     key_max3 = key_list[val_list.index(sorted_val_list[-3])]
    # else:
    #     key_max3 = 0

    if key_max1 in lst[it].tab:    
        v0 = lst[it].tab[key_max1]
    else:
        v0 = 0.

    if key_max2 in lst[it].tab:    
        v1 = lst[it].tab[key_max2]
    else:
        v1 = 0.
    v2 = 0 # lst[it].tab[maximum3]

    R=1
    G=1
    B=1

    colors=[[None for i in range(53) ]]
    BitLabel=[[str(i) for i in range(53) ]]
    V=v0+v1+v2
    alpha=1.
    for i in range(53):
        if (v0/V>1):
            red = R
            R -= 0.018
        else:
            red = 0

        if (v1/V>1):
            green = G
            G -= 0.018
        else:
            green = 0

        if (v2/V>1):
            blue = B
            B -= 0.018
        else:
            blue = 0

        V /= 2
        colors[0][i] = tuple((red,green,blue,alpha)) 

    V=v0+v1+v2
    for i in range(52,0,-1):
        if (V>1):
            colors[0][i]=tuple((colors[0][i][0],colors[0][i][1],colors[0][i][2], 0.4))
        V/=2
    
 
    colorlegend1 = mpatches.Patch(color='red', label=str(key_max1))
    colorlegend2 = mpatches.Patch(color='green', label=str(key_max2))
  #  colorlegend3 = mpatches.Patch(color='blue', label=str(key_max3))
    colorlegend4 = mpatches.Patch(color='white', label="Cancelled bit")

    fig=plt.figure(linewidth=2,
                #edgecolor='steelblue', # 
                #facecolor='skyblue',
                tight_layout={'pad':1},
                figsize=(15,3))

    the_table = plt.table(cellText=BitLabel, cellColours=colors, loc='center', cellLoc='center')
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(10)


    #the_table.patch.set_facecolor('blue')

    #the_table = ax.table(cellColours=colors,loc='center')
    #the_table = plt.table(rowColours=([(0.1,0.1,0.8),(0.8,0.1,0.1)]),cellText=[['toto'],['titi']], loc='center')
    plt.box(on=None)
    plt.legend(handles=[colorlegend1, colorlegend2, colorlegend4])
    plt.axis('tight')
    plt.axis('off')
    plt.suptitle("Evolution of T1 weight in the computed results of PI with Archimede's iterations")
    plt.title("Iteration N°: "+str(it)+" Value: "+str(lst[it].val))

    #plt.draw()
    plt.savefig("pi_archimedes"+str(it)+".png")
    plt.close(fig)
    #plt.show()
    #im = plt.imshow(animated=True)
    #ims.append([im])

# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
# ani.save("pi_archimedes.gif", writer="imagemagick")
plt.show()




