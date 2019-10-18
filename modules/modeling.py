#------------------------------------------------------------------------------------
# Title: Model Codes
# Author: Bruno Lima de Freitas
# Description: Source codes to model objects.
#------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------------
def multi_prism( index, n_prism, xmin, xmax, zmin, zmax):
    '''
    This function generates a group of juxtaposed prisms in a limited space and plot 
    some or all of them. All the generated prisms presents the same horizontal and 
    vertical extension. Futhermore, all the prisms presents the same top level and 
    bottom level.
    
    Inputs:
    index : (list) ==> indexes of plotted prisms    
    n_prism : (interger) ==> number of generated prisms
    xmin : (interger/float) ==> minimum value of the horizontal limits of the space
    where prisms will be plotted
    xmax : (interger/float) ==> maximum value of the horizontal limits of the space
    where prisms will be plotted
    zmin : (interger/float) ==> minimum value of the vertical limits of the space
    where prisms will be plotted
    zmax : (interger/float) ==> maximum value of the vertical limits of the space
    where prisms will be plotted
    
    '''
    # Confirming conditions for prisms's selection
    if max(index) >= n_prism:
        raise ValueError("Remove incorrect indexes from index list")
    
    # Defining horizontal coordinates of prisms's center
    x = np.linspace(xmin, xmax, n_prism)             
    
    # Defining horizontal and vertical coordinates of prisms's vertex
    xmed = ( x[0] + x[1] ) / 2.0
    x_prism = abs( x[0] - ( xmed ) ) # distance of prism's center to one of your horizontal limits     
    
    x_coord = [] # horizontal coordinates of each vertex of each prism
    z_coord = [] # vertical coordinates of each vertex of each prism
    # x_coord[i] and z_coord[i] corresponds the coordinates of prism which index is i
    for i in range (n_prism):
        x_coord.append([x[i] - x_prism, x[i] + x_prism, x[i] + x_prism, \
                        x[i] - x_prism, x[i] - x_prism])
   
        z_coord.append([zmin, zmin, zmax, zmax, zmin])            
    
    # Plotting the selected prisms
    figure, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(17, 4), facecolor='w')
                
    for i in index:
        ax1.plot(x_coord[i], z_coord[i], ".-r", linewidth=2)

        ax1.fill_between(x_coord[i], z_coord[i], facecolor='red', alpha=0.1)

    ax1.set_ylabel('Depth (m)', fontsize=15)
    ax1.set_xlabel('Coordinate x (m)', fontsize=15)
    ax1.invert_yaxis()
    ax1.label_outer()
    ax1.set_xlim(xmin - x_prism, xmax + x_prism)
    plt.show( )

#####################################################################################