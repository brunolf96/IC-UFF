#------------------------------------------------------------------------------------
# Title: Model Codes
# Author: Bruno Lima de Freitas
# Description: Source codes to model objects.
#------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
#------------------------------------------------------------------------------------
def multi_prism(n_prism, xmin, xmax, zmin, zmax, index='all'):
    '''
    This function generates a group of juxtaposed prisms in a limited space and plot some or all of them. All the generated prisms presents the same horizontal extension. Futhermore, all the prisms, not necessarily, presents the same the top level and bottom level.
    
    Inputs:
    n_prism : (interger) ==> Number of generated prisms
    xmin : (interger/float) ==> Minimum value of the horizontal limits of the space where prisms will be plotted
    xmax : (interger/float) ==> Maximum value of the horizontal limits of the space where prisms will be plotted
    index : (list/string) ==> Default value is 'all' and represents that all the prisms genarated between xmin and xmax will be used. If only some of the defined prism will be used, index should be a list that contains the indexes of those prisms 
    zmin : (1D array) ==> Minimum values of the vertical limits of each prism generated that will be plotted. If z_const is 'yes', zmin array should contain just one value
    zmax : (1D array) ==> Maximum values of the vertical limits of each prism generated that will be plotted. If z_const is 'yes', zmax array should contain just one value
    
    '''       
    # Confirming types of variables
    if type(n_prism) != int:
        raise ValueError('Variable n_prism should be a interger')
    if type(xmin) != int and type(xmin) != float:
        raise ValueError('Variable xmin should be a interger or a float')
    if type(xmax) != int and type(xmin) != float:
        raise ValueError('Variable xmax should be a interger or a float')
    if isinstance(zmin, np.ndarray) == False:
        raise ValueError('Variable zmin should be a array')
    if isinstance(zmax, np.ndarray) == False:
        raise ValueError('Variable zmax should be a array')
    if index != 'all' and type(index) != list:
        raise ValueError('Variable index should be a string with value ''all'' or should be a list') 
    
    # Confirming conditions for prisms's selection
    if index != 'all':
        if max(index) >= n_prism:
            raise ValueError('Remove incorrect indexes from index list')
        if min(index) < 0:
            raise ValueError('Remove incorrect indexes from index list')
        
    # Confirming conditions for zmax and zmin:
    if len(zmax) > 1 and len(zmin) > 1:
        if index == 'all':
            if len(zmax) != n_prism:
                raise ValueError('Variable zmax presents a number of values different from the number of selected prisms')            
            if len(zmin) != n_prism:
                raise ValueError('Variable zmin presents a number of values different from the number of selected prisms')   
        else:
            if len(zmax) != len(index):
                raise ValueError('Variable zmax presents a number of values different from the number of selected prisms')            
            if len(zmin) != len(index):
                raise ValueError('Variable zmin presents a number of values different from the number of selected prisms') 
    else:
        if len(zmax) == 1:
            if len(zmin) != 1:
                raise ValueError('Verify variables zmin and zmax')
        if len(zmin) == 1:
            if len(zmax) != 1:
                raise ValueError('Verify variables zmin and zmax')

    # Defining the indexes of prisms that will be used
    if index == 'all':
        index = []
        for i in range (n_prism):
            index.append(i)
    
    # Defining horizontal coordinates of prisms's center
    x = np.linspace(xmin, xmax, n_prism)             
    
    # Defining horizontal and vertical coordinates of prisms's vertex
    xmed = ( x[0] + x[1] ) / 2.0
    x_prism = abs( x[0] - ( xmed ) ) # distance of prism's center to one of your horizontal limits     
    
    x_coord = [] # horizontal coordinates of each vertex of each prism
    z_coord = [] # vertical coordinates of each vertex of each prism
    # x_coord[i] and z_coord[i] corresponds the coordinates of prism which index is i
    aux = 0 # variable created to permit the correct use of variables zmin and zmax
    for i in range (n_prism):
        x_coord.append([x[i] - x_prism, x[i] + x_prism, x[i] + x_prism, \
                        x[i] - x_prism, x[i] - x_prism])
   
        if len(zmax) == 1 and len(zmin) == 1:
            z_coord.append([zmin[aux], zmin[aux], zmax[aux], zmax[aux], zmin[aux]])
        else:
            if i in index:
                z_coord.append([zmin[aux], zmin[aux], zmax[aux], zmax[aux], zmin[aux]])
                aux += 1
            else:
                z_coord.append([0, 0, 0, 0, 0])
    
    # Plotting the selected prisms
    figure, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(17, 4), facecolor='w')
                
    for i in index:
        ax1.plot(x_coord[i], z_coord[i], ".-k", linewidth=2)

    ax1.set_ylabel('Depth (m)', fontsize=15)
    ax1.set_xlabel('Coordinate x (m)', fontsize=15)
    ax1.invert_yaxis()
    ax1.label_outer()
    ax1.set_xlim(xmin - x_prism, xmax + x_prism)
    plt.show( )
    
    return x_prism, x_coord, z_coord

#####################################################################################