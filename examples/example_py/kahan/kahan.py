# Importing NumPy Library
import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
import seaborn as sns               # For the automatic handling of gradient for the heatmap  
import sys
sys.path.append('../../..')

from random import *
from chnbr import ChNbr

N=20

# 1. IEEE-754
x = 1.60631924

def cf(x):
    V = 4. - (3.*(x-2.)*((x-5.)*(x-5.)+4.)) / (x + (x-2.0)*(x-2.0)*((x-5.)*(x-5.) + 3.) )
    return V

def rp(x):
    V = (622. - x*(751. - x*(324. - x*(59. - 4.*x)))) / (112.-x*(151.-x*(72.-x*(14.-x))))
    return V

for i in range(1,N):
    dx = i*2**(-53)
    kh = cf(x+dx)-rp(x)
#    print(i, dx, kh)


# 2. CHNBR
# lst : List of chromatic numbers
lst=[]
x = ChNbr(1.60631924, follow=True)
eps = ChNbr(2**(-53), follow=True)  
c2= ChNbr(2., follow=True) 
c3= ChNbr(3., follow=True) 
c4= ChNbr(4., follow=True) 
c5= ChNbr(5., follow=True) 
c14= ChNbr(14., follow=True)
c59= ChNbr(59., follow=True)
c72= ChNbr(72., follow=True)
c112= ChNbr(112., follow=True)
c151= ChNbr(151., follow=True)
c324= ChNbr(324., follow=True)
c622= ChNbr(622., follow=True)
c751= ChNbr(751., follow=True)

def cf(x):
    V = c4 - (c3*(x-c2)*((x-c5)*(x-c5)+c4)) / (x + (x-c2)*(x-c2)*((x-c5)*(x-c5) + c3) )
    return V

def rp(x):
    V = (c622 - x*(c751 - x*(c324 - x*(c59 - c4*x)))) / (c112-x*(c151-x*(c72-x*(c14-x))))
    return V

for i in range(1,N):
    dx = i*eps
    kh = cf(x+dx)-rp(x)
    lst.append(kh)




# plot results
for it in range(0,N-1):
    # list out keys and values separately
    key_list = list(lst[it].tab.keys())
    val_list = list(lst[it].tab.values())

    sorted_val_list = sorted(val_list)

    key_max1 = key_list[val_list.index(sorted_val_list[-1])]
    if len(val_list)>1:
        key_max2 = key_list[val_list.index(sorted_val_list[-2])]
    else:
        key_max2 = 0

    # 1st largest => RED
    if key_max1 in lst[it].tab:    
        v0 = lst[it].tab[key_max1]
    else:
        v0 = 0.

    # 2nd largest => GREEN
    if key_max2 in lst[it].tab:    
        v1 = lst[it].tab[key_max2]
    else:
        v1 = 0.

    # The rest => BLUE
    v2=0

    for k, v in lst[it].tab.items():
        if k!=key_max1 and k!=key_max2:
            v2 += v

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
        colors[0][i] = tuple((red, green, blue, alpha)) 

    V=v0+v1+v2
    for i in range(52,-1,-1):
        if (V>1):
            colors[0][i]=tuple((colors[0][i][0],colors[0][i][1],colors[0][i][2], 0.4))
        V/=2
    
 
    colorlegend1 = mpatches.Patch(color='red', label=str(key_max1))
    colorlegend2 = mpatches.Patch(color='green', label=str(key_max2))
    colorlegend3 = mpatches.Patch(color='blue', label="Other")
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
    plt.legend(handles=[colorlegend1, colorlegend2, colorlegend3, colorlegend4])
    plt.axis('tight')
    plt.axis('off')
    plt.suptitle("Kahan's Iteration")
    plt.title("Iteration N°: "+str(it)+" Value: "+str(lst[it].val))

    #plt.draw()
    plt.savefig("Kahan"+str(it)+".png")
    plt.close(fig)
    #plt.show()
    #im = plt.imshow(animated=True)
    #ims.append([im])

# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
# ani.save("pi_archimedes.gif", writer="imagemagick")
plt.show()
