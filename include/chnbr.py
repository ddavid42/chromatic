""" 
    Module To follow relative weight of variables using chromatic analysis
    Define ChNbr type and associated methods 
    Example:
    
    c1=ChNbr(1, follow=True)
    c3=c1+2
    c3 -= 2
    c3 += 10

    c4=c1+3
    print(c1)
    print(c3)
    print(c4)
"""
from math import *
from re import M
import copy
import numpy as np

# This version correspond to an online computation of an approximation of the Condition Number

# Chromatic number class
# Global variables
""" 
Idx : 
    Global index used to store value to track 
ForgiveFactor : 
    Factor between the largest and the lowest number used to no longer track specifically
    variables V_i those relative weight is such that 
    |V_i|/|V_m| < ForgiveFactor with V_m corresponding the largest variable
"""
Idx = 0      
ForgiveFactor = 1e-10

def get_newindex():
    """
    Global function returning a unique global identifier
    """
    global Idx
    Idx += 1 
    return Idx

class ChNbr:
    """
    Class used to represent chromatic numbers. 
    Chromatic numbers is a way to track the relative weight of variables trough computations.
    """
    val = 0.0
    tab = {}

    def __init__(self, val=0.0, tab={-1:1.0}, follow=False):
        """
        Create a new chromatic number with the following argument
        val: input value
        tab: dictionnary storing the impact of each followed value
        follow = define if this value has to be followed
        """
        
        self.val = val
        if (follow):
            idx = get_newindex()
            self.tab = {idx:1.0}
            self.idx = idx
        else:
            self.idx = -1
            self.tab = dict()
            max_key = max(tab.values())
            for key in tab:
                if (tab[key]/max_key) > ForgiveFactor :
                    self.tab[key] = tab[key]
                else:
                    if -2 in self.tab:
                        self.tab[-2] += tab[key]
                    else:
                        self.tab[-2]  = tab[key]


    def __repr__(self):
        s = "value: " + str(self.val) + ", "
        s+= "Idx  : " + str(self.idx)
        for k,v in self.tab.items():
            s+= ", [" + str(k) + "," + str(v) + "]"
        return s

    ## Operations
    def __str__(self):
        """ Convert a Chromatic Number to a string """
        s = "value: " + str(self.val) + "\n"
        s+= "Idx  : " + str(self.idx) + "\n"
        for k,v in self.tab.items():
            s+= "[" + str(k) + "," + str(v) + "]\n"
        return s

    def __neg__(self):
        """ -(A) with A a chromatic number """
        return ChNbr(-(self.val), self.tab)

    def __add__(self, other):
        """ (A)+(B) with A a chromatic number """

        lhs = self
        
        if type(other) is not ChNbr:
            if not np.isscalar(other):
                return other+self
            rhs = ChNbr(other)
        else:
            rhs = other

        # Value
        left  = lhs.val
        right = rhs.val
        res = left + right

        # New empty Dictionary
        tab=dict()
        
        # WARNING: when res=0 we set coef=1 to avoid division by 0 arbitrarily but 
        # infty could be a more logical option for 'left' and 'right' member
        if res != 0:
            coef = fabs(left/res)
        else:
            coef = float('inf')
        for i,t in lhs.tab.items():
            tab[i] = coef*t

        if res != 0:
            coef = fabs(right/res)
        else:
            coef = float('inf')  

        for i,t in rhs.tab.items():
            if i in tab:
                tab[i] += coef * t
            else:
                tab[i] = coef * t

        # Create the number
        nbr = ChNbr(res,tab)
        return nbr

    def __radd__(self,other):
        """ (A)+(B) with B a chromatic numbers and A is not """
        return self + other

    def __sub__(self, other):
        """ (A)-(B) with A a chromatic number """
        return self + (-other)

    def __rsub__(self, other):
        """ (A)-(B) with B a chromatic numbers and A is not """
        return (-self) + other

    def __mul__(self, other):
        """ (A)*(B) with A a chromatic number """
     
        lhs = self
        
        if type(other) is not ChNbr:
            # Handle special cases
            if not np.isscalar(other):
                return other*self
            elif other==0:
                return 0
            elif other==1:
                return ChNbr(self.val, self.tab)
            else:
                rhs = ChNbr(other)
        else:
            rhs = other

        # Value
        res = lhs.val * rhs.val

        # Dictionary
        tab=dict()
        for i,t in rhs.tab.items():
            tab[i]=t/2.0

        for i,t in lhs.tab.items():
            if i in tab:
                tab[i]+=t/2.0
            else:
                tab[i]=t/2.0

        # Create the number
        nbr = ChNbr(res,tab)
        return nbr

    def __rmul__(self, other):
        """ (A)*(B) with B a chromatic numbers and A is not """
        return self*other        

    def __truediv__(self, other):
        """ (A)/(B) with A a chromatic number """
 
        lhs = self
        
        if type(other) is not ChNbr:
            # Handle special cases
            if other==0:
                return float('inf')
            elif other==1:
                return ChNbr(self.val, self.tab)
            else:
                rhs = ChNbr(other)
        else:
            rhs = other

        # Value
        res = lhs.val / rhs.val

        # Dictionary
        tab=dict()
        for i,t in rhs.tab.items():
            tab[i]=t/2.0

        for i,t in lhs.tab.items():
            if i in tab:
                tab[i]+=t/2.0
            else:
                tab[i]=t/2.0

        # Create the number
        nbr = ChNbr(res,tab)
        return nbr

    def __rtruediv__(self, other):
        """ (A)/(B) with B a chromatic numbers and A is not """
        
        if type(self) is not ChNbr:
            # Handle special cases
            if self==0:
                return (0/other.val)

        return ((ChNbr(1)/self) * other) 



    # Comparisons
    def __eq__(self, other):
        """(A)==(B)"""
        if type(other) is ChNbr:
            return (self.val == other.val) 
        else:
            return (self.val == other)

    def __ne__(self, other):
        """(A)!=(B)"""
        if type(other) is ChNbr:
            return (self.val != other.val) 
        else:
            return (self.val != other)

    def __lt__(self, other):
        """(A) < (B)"""
        if type(other) is ChNbr:
            return (self.val < other.val) 
        else:
            return (self.val < other)

    def __le__(self, other):
        """(A) <= (B)"""
        if type(other) is ChNbr:
            return (self.val <= other.val) 
        else:
            return (self.val <= other)

    def __gt__(self, other):
        """(A) > (B)"""        
        if type(other) is ChNbr:
            return (self.val > other.val) 
        else:
            return (self.val > other)

    def __ge__(self, other):
        """(A) >= (B)"""
        if type(other) is ChNbr:
            return (self.val >= other.val) 
        else:
            return (self.val >= other)

    @classmethod
    def sqrt(nbr1,nbr2):
        nbr = ChNbr(sqrt(nbr2.val))
        nbr.tab = nbr2.tab.copy()
        return nbr 

    def tanh(self):
        own_contrib = abs(0.5*self.val*(1.0 - tanh(self.val)**2))
        nbr = ChNbr(tanh(self.val))
        nbr.tab = self.tab.copy()
        if own_contrib < abs(nbr.val):
            for i, t in nbr.tab.items():
                nbr.tab[i] *= (own_contrib)/abs(nbr.val)
            if not(-1 in nbr.tab):
                nbr.tab[-1] = 0.0
            nbr.tab[-1] += (abs(nbr.val) - own_contrib)/abs(nbr.val)
        return nbr 

    def __pow__(self, exponent : int):
        nbr = ChNbr(pow(self.val, exponent))
        nbr.tab = self.tab.copy()
    
    def __pow__(self, exponent : int):
        nbr = ChNbr(pow(self.val, exponent))
        nbr.tab = self.tab.copy()
        return nbr 