# -----------------------------------------------------------------------------------------------------------------------------------------
# Title: Gravity Codes
# Description: Source codes for calculate gravity.
# -----------------------------------------------------------------------------------------------------------------------------------------

# Import Python libraries:
import numpy as np 

###########################################################################################################################################
def g_sphere(x, z, sphere1, component='z'):
    '''    
    This function calculates all components of gravity attraction produced by a solid point mass and returns the one associated to the required one.
    This is a Python implementation for the subroutine presented in Blakely (1995). On this function, there are received the value of the initial
    and final observation points (X and Y) and the properties of the sphere.
       
    Inputs:
    x - numpy array - observations in x directions (meters)
    y - numpy array - observations in y directions (meters)
    component - string - the required component to be calculated
    sphere1 - list - elements of the sphere: [x_center(meters), z_center(meters), mass(kg)]
    
    Output:
    g - numpy array - the required component for the gravity in mGal. Size of gz is the same as x and z observations    
    '''
    
    # Stablishing some conditions
    if x.shape != z.shape:
        raise ValueError("All inputs must have same shape!")
    
    # Definition for some constants
    G = 6.673e-11 # SI
    si2mGal = 100000.0
    
    
    # building a "transposed" list for correct usage of enumerate loop:
 
    sphere=[]
    for i in range( len(sphere1[0]) ): # sphere1[0] contains of all masses considered.
        sphere.append([ sphere1[0][i], sphere1[1][i], sphere1[2][i] ])

    # Setting the initial value for gravity:
    g = 0.
    gg = 0.0
    xs = np.zeros( len(sphere) )
    zs = np.zeros( len(sphere) )
    ms = np.zeros(len(sphere) )
    # loop for all point masses in list sphere
    for i,j in enumerate(sphere):
        #print (i,j)
        xs[i] = j[0]
        zs[i] = j[1]
        ms[i] = j[2]
    
        # Setting position-vector: 
        dx = xs[i] - x
        dz = zs[i] - z
    
        # properties of the sphere:
        mass = ms[i]
     
        # Compute the distance
        r = np.sqrt(dx**2 + dz**2)
    
        if component=='z':
            # Compute the vertical component 
            g = mass * dz / (r**3)
            g *= G*si2mGal
        elif component =='x':
            # Compute the vertical component 
            g = mass * dx / (r**3)
            g *= G*si2mGal
        
        gg += g
        # Return the final outpu
    return gg

###########################################################################################################################################