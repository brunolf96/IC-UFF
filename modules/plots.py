# -----------------------------------------------------------------------------------------------------------------------------------------
# Title: Grav-Mag Codes
# Author: Bruno Lima 
# Colaborador: Rodrigo Bijani 
# Description: Source codes for plotting images.
# -----------------------------------------------------------------------------------------------------------------------------------------

# Import Python libraries:
import numpy as np
import matplotlib.pyplot as plt
import pylab as py

# bibliotecas para pintar prismas!
from matplotlib.path import Path
from matplotlib.patches import PathPatch

###########################################################################################################################################
def plots_rectangles(xo, zo, p, color='black', paint=False, var=None, name=None, repeat=False):
    """
    Plot a group of rectangles with the same width in a juxtaposed framework following a number of observations.

    Inputs:
    
    * xo : (1D array) ==> x coordinate of observations;
    * zo : (1D array) ==> z coordinate of observations;   
    *  p : (1D array) ==> difference between top and bottom of a specific rectangle;     
    * color : (string) ==> color of the rectangles; 
    
    * paint : (string) ==> If paint == 'True', all the rectangles will be painted based on the values of the var. If paint = 'False', the rectangles won't be painted;
    
    * var : (list) ==> If you paint the rectangles, you should pay attention on this variable. Values related to a property or physical greatness of the layers on the subsurface. These values have to be given in a specific form that will depend if repeat == 'True' or repeat == 'False'. 
    The examples given below should be considered if repeat == 'False':
       - If you have just one rectangle to paint: var = [ [13,1,10,31,54] ]
       - If you have two rectangles to paint: var = [ [14,25,51,32,61], [2,4,6,7,1] ]
       - If you have three rectangles to paint: var = [ [14,16,71,10,11], [21,22,52,62,29], [32,13,43,23,39] ]
    The examples given below should be considered if repeat == 'True':
       - If you have just one rectangle to paint: var = [ [13] ]
       - If you have two rectangles to paint: var = [ [14], [2] ]
       - If you have three rectangles to paint: var = [ [71], [21], [32] ];
    The numbers of rectangles to be paint will define the numbers of lists inside the list;
    
    * name : (string) ==> If you paint the rectangles, you should pay attention on this variable. It refers to the name of the property or physical greatness of the layers on the subsurface that is in study. If the rectangles will be painted, you can provide the name of this property or physical greatness as you can see in the examples below. This name will be show on the colorbar.
    For example:
        - name = 'Depth $(m)$'
        - name = 'Density constrast $(g/cm^3)$';
    
    * repeat : (string) ==> If you paint the rectangles, you should pay attention on this variable. If repeat == 'True', it means the values of var are constant in depth and just one value have to be given for each rectangle. If repeat == 'False', it means the values of var changes in depth and five values of var have to be given for each rectangle. In the case of repeat == 'False', each value of var for a specific rectangle are related with each point that forms the rectangle. The form that var will be given to this function will depend of repeat and it is clearly detailed on the description of var;
    
    Outputs:
    
    * x_plot : (list) ==> It is a list that contains other lists where each of them provides the x coordinate of the points that form a specific rectangle that has been plotted;
    * z_plot : (list) ==> It is a list that contains other lists where each of them provides the z coordinate of the points that form a specific rectangle that has been plotted;

    """
    # metade da espessura horizontal do primeiro retangulo em funcao da posição de observacao (1 prisma por observacao):
    # localizacao do primeiro prisma em funcao da posicao de observacao (1 prisma por observacao):
    xmed = ( xo[0] + xo[1] ) / 2.0
    x_prisma = abs( xo[0] - ( xmed ) ) # esse valor será usado durante a etapa de visualização gráfica para a possibilizar a geração dos retangulos
    
    # criando listas que armazenarão as coordenadas usadas nos gráficos produzidos: 
    nobs = len(xo)
    x_plot = []
    z_plot = []
    for i in range(nobs):
        x_plot.append( [ xo[i] - x_prisma, xo[i] + x_prisma, xo[i] + x_prisma, xo[i] - x_prisma, xo[i] - x_prisma ] )
        # O valor 0 foi usado como referência
        # valores acima ou abaixo de 0 para z0, simbolizando um excesso ou falta de relevo, respectivamente, em relação a referência
        if zo[i] < 0:
            z_plot.append( [ zo[i], zo[i], 0 + p[i] , 0 + p[i], zo[i] ] )
        else:
            z_plot.append( [ zo[i], zo[i], zo[i] + p[i] , zo[i] + p[i], zo[i] ] )
    
    if paint == False:
        
        # Visualização gráfica:
        plt.figure( figsize=(10,10) )
        plt.plot(xo,zo,'vr')
    
        for i in range(nobs):
            plt.plot( x_plot[i], z_plot[i], color) # visualização dos retângulos
    
        plt.grid()
        plt.xlim( [ xo[0] - 2 * x_prisma , xo[nobs - 1] + 2 * x_prisma ] ) # retirar esse comando, interfere na visualização
        plt.gca().invert_yaxis()
        plt.show()
    
    else:
        xo2 = x_plot
        zo2 = z_plot
        p2 = p
        var2 = var
        name2 = name
        repeat2 = repeat
        plots_paint_rectangles(xo2, zo2, p2, var=var2, name=name2, repeat=repeat2)
        
    return x_plot, z_plot

###########################################################################################################################################
def plots_paint_rectangles(xo, zo, p, var=None, name=None, repeat=False):
    """
    Paint a group of rectangles with the same width in a juxtaposed framework.;5

    Inputs:
    
    * xo : (list) ==> It is a list that contains other lists where each of them provides the x coordinate of the points that form a specific rectangle that has been plotted;
    * zo : (list) ==> It is a list that contains other lists where each of them provides the z coordinate of the points that form a specific rectangle that has been plotted;
    *  p : (1D array) ==> difference between top and bottom of a specific rectangle;
    * var : (list) ==> Values related to a property or physical greatness of the layers on the subsurface. These values have to be given in a specific form that will depend if repeat == 'True' or repeat == 'False'. 
    The examples given below should be considered if repeat == 'False':
       - If you have just one rectangle to paint: var = [ [13,1,10,31,54] ]
       - If you have two rectangles to paint: var = [ [14,25,51,32,61], [2,4,6,7,1] ]
       - If you have three rectangles to paint: var = [ [14,16,71,10,11], [21,22,52,62,29], [32,13,43,23,39] ]
    The examples given below should be considered if repeat == 'True':
       - If you have just one rectangle to paint: var = [ [13] ]
       - If you have two rectangles to paint: var = [ [14], [2] ]
       - If you have three rectangles to paint: var = [ [71], [21], [32] ];
    The numbers of rectangles to be paint will define the numbers of lists inside the list;
    
    * name : (string) ==> It refers to the name of the property or physical greatness of the layers on the subsurface that is in study. If the rectangles will be painted, you can provide the name of this property or physical greatness as you can see in the examples below. This name will be show on the colorbar.
    For example:
        - name = 'Depth $(m)$'
        - name = 'Density constrast $(g/cm^3)$';
    
    * repeat : (string) ==> If repeat == 'True', it means the values of var are constant in depth and just one value have to be given for each rectangle. If repeat == 'False', it means the values of var changes in depth and five values of var have to be given for each rectangle. In the case of repeat == 'False', each value of var for a specific rectangle are related with each point that forms the rectangle. The form that var will be given to this function will depend of repeat and it is clearly detailed on the description of var;
    
    Outputs:

    """
    plt.figure( figsize=(10,10) )
    
    # metade da espessura horizontal do primeiro retangulo em funcao da posição de observacao (1 prisma por observacao):
    # localizacao do primeiro prisma em funcao da posicao de observacao (1 prisma por observacao):
    xmed = ( xo[0][0] + xo[1][0] ) / 2.0
    x_prisma = abs( xo[0][0] - ( xmed ) ) # esse valor será usado durante a etapa de visualização gráfica para a possibilizar a geração dos retangulos
    
    nobs = len(xo)
    if repeat == True:
        for i in range ( nobs ):
            for j in range ( len(xo[0]) - 1): # utilizamos xo[0] porque precisamos associar um valor de var para cada ponto que usaremos para plotar o retângulo mas qualquer outro índice no lugar de zero é válido
                var[i].append( var[i][0] )
                
    for i in range ( nobs ):
        for j in range( len(xo[0]) - 1 ): # utilizamos xo[0] porque precisamos associar um valor de var para cada ponto que usaremos para plotar o retângulo mas qualquer outro índice no lugar de zero é válido
            value = var[i][j]
            if j == 0 and i == 0:
                var_min = value
                var_max = value
            else:
                if value > var_max:
                    var_max = value
                if value < var_min:
                    var_min = value       
        print('i =',i,'=>',var[i]) # apenas usado para confirmar mais facilmente a escala de cor até ela ser ajustada corretamente
       
    for i in range ( nobs):
        if zo[i][0] < 0:
            plt.plot( xo[i], zo[i], 'k.-')
            zp = np.array(zo[i])
            if i == 0:
                zmin = zo[i][0] 
                zmax = 0 + p[i]
            else:
                if 0 + p[i]  > zmax:
                    zmax = 0 + p[i]
                if zo[i][0]  < zmin:
                    zmin = zo[i][0] 
        else:
            plt.plot( xo[i], zo[i], 'k.-')
            zp = np.array(zo[i])
            if i == 0:
                zmin = zo[i][0] 
                zmax = zo[i][0] + p[i]
            else:
                if zo[i][0] + p[i]  > zmax:
                    zmax = zo[i][0] + p[i]
                if zo[i][0]  < zmin2:
                    zmin = zo[i][0] 
    
        xp = np.array( xo[i] )    
    
        path = Path(np.array([xp,zp]).T)
        patch = PathPatch(path, facecolor='none')
    
        plt.gca().add_patch(patch)
        fs = 18 # font size for the label
        
        var_part = np.array( var[i] ) # Isolando apenas os dados de var que serão usados nesse for
        
        im = plt.imshow(var_part.reshape(np.size(zp),1), cmap='RdBu',interpolation="bicubic", 
                    vmin=var_min, vmax=var_max,
                    origin='lower',extent=[min(xp), max(xp), min(zp), max(zp)],
                    aspect="auto", clip_path=patch, clip_on=True)
   
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(name, fontsize=fs)
    
    plt.xlim( [ xo[0][0] - 2 * x_prisma , xo[nobs - 1][1] + 2 * x_prisma ] ) # retirar esse comando, interfere na visualização
    plt.ylim(zmax + 2, zmin - 2) # retirar esse comando, interfere na visualização

    plt.grid()
    plt.show()
        
###########################################################################################################################################
        
        
        
        
        
        