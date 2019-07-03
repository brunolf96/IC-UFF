# -----------------------------------------------------------------------------------------------------------------------------------------
# Title: Gravity Codes
# Description: Source codes for calculate gravity.
# -----------------------------------------------------------------------------------------------------------------------------------------

# Import Python libraries:
import numpy as np 

###########################################################################################################################################
def g_prism(x,z,prism):    
    '''
    This function calculates the vertical component of gravity attraction produced by a prism 
    
    (Telford, 1981)
    Inputs: x,z = arrays with cartesian coordinates in meters;
    rod: list with the following elements in this specific order: prism[x_prism, h1_prism, L, A, rho ]:
    x_prism      = horizontal distance from origin (meters)
    h1_prism      = vertical distance from origin to top (meters)
    L         = length of the thin prism (meters) (never zero!)
    rho       = density of the prism (kg/m3)
    A         = cross section of the prism (m2)
    component = the gravitational component to be computed (only z component is available yet) 
    
    Output:
    g - numpy array - the required component for the gravity in mGal. Size of gz is the same as x and z observations    
    '''

    #saving the inputs of prism:
    x_prism  = prism[0]
    h1_prism = prism[1]
    L_prism  = prism[2]
    h2_prism  = L_prism - h1_prism
    A        = prism[3]
    rho      = prism[4]
    
    #seting the variables to make the calclation?
    dx = x - x_prism
    dy = np.zeros(len(x))
    h1 = z - h1_prism
    h2 = z - h2_prism
    
    #setting some constans:
    G = 6.673e-11
    si2mGal = 100000.0
    
    #making the comptation of g
    term0 = (G*rho*A)
    
    term11 = ( dx**2 + dy**2 + h1**2)
    term1 = 1/(term11**(0.5))
    
    term21 = ( dx**2 + dy**2 + h2**2)
    term2 = 1/(term21**(0.5))

    g = term0 * abs(term1 - term2)
    
    return g*si2mGal

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
def g_prism2D_Bruno(x,z,prism):
    '''
    This function calculates the vertical component of gravity attraction produced by a rectangle
    
    Inputs: 
    x = horizontal cartesian coordinate in meters of the observation;
    z = vertical cartesian coordinate in meters of the observation;
    rod: list with the following elements in this specific order: prism[x_prism, x_med, z_prism, L, rho ]:
    x_prism = horizontal cartesian coordinate of the center of the rectangle (meters)
    z_prism = vertical cartesian coordinate of the center of the rectangle (meters)
    x_med   = half of horizontal length of rectangle (meters)
    L = vertical length of the rectangle (meters)
    rho       = density of the prism (kg/m3)
    component = the gravitational component to be computed (only z component is available yet) 
    
    Output:
    g - numpy array - the required component for the gravity in mGal. Size of gz is the same as x and z observations    
    '''
    
    # saving the inputs of prism:
    x_prism = prism[0]
    x_med = prism[1]
    z_prism = prism[2]
    L = prism[3]
    rho = prism[4]
    
    # constructing the lengths of the rectangle
    x_length = 2 * x_med
    y_length = 1
    z_length = L
    
    # calculating the volume of the rectangle
    volume = x_length * y_length * z_length
    
    # distancia entre as observacoes e a fonte:
    dx = x_prism - x
    dz = z_prism - z
    r = np.sqrt(dx**2 + dz**2)
    
    # essencial informations
    gamma = 6.67e-11 # no SI (m3 / s2 kg)
    si2mGal = 100000.0 # converter para mGal (unidade convencional)

    # calculating the gravimetric contribution of the rectangle
    g = (gamma * rho * volume * dz ) / r**3
    g *= si2mGal
    
    # Return the final output
    return g

###########################################################################################################################################