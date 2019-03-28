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
def plots_rectangles(xo, zo, p, ref=0, color='black'):
    """
    Plot a group of rectangles with the same width in a juxtaposed framework following a number of observations.

    Inputs:
    
    * xo : (1D array) ==> x coordinate of observations;
    * zo : (1D array) ==> z coordinate of observations;   
    *  p : (1D array) ==> difference between top and bottom of a specific rectangle;     
    * color : (string) ==> color of the rectangles; 
    * ref : (int) ==> represents the value of reference level that will be used to establish the position of the top of the rectangles;
    
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
            z_plot.append( [ zo[i], zo[i], ref + p[i] , ref + p[i], zo[i] ] )
        else:
            z_plot.append( [ zo[i], zo[i], zo[i] + p[i] , zo[i] + p[i], zo[i] ] )
        
    # Visualização gráfica:
    plt.figure( figsize=(10,10) )
    plt.plot(xo,zo,'vr')
    
    for i in range(nobs):
        plt.plot( x_plot[i], z_plot[i], color) # visualização dos retângulos
    
    plt.grid()
    plt.xlim( [ xo[0] - 2 * x_prisma , xo[nobs - 1] + 2 * x_prisma ] ) # retirar esse comando, interfere na visualização
    plt.gca().invert_yaxis()
    plt.show()
        
    return x_plot, z_plot

###########################################################################################################################################
def plots_paint_rectangles(xo, zo, p, ref=0, n_var=None, var=None, name=None, cmap='RdBu_r'):
    """
    Paint a group of rectangles with the same width in a juxtaposed framework.

    Inputs:
    
    * xo : (1D array) ==> x coordinate of observations;
    * zo : (1D array) ==> z coordinate of observations;       
    *  p : (1D array) ==> difference between top and bottom of a specific rectangle; 
    * cmap : (string) ==> colormap instance or registered colormap name;
    * ref : (int) ==> represents the value of reference level that will be used to establish the position of the top of the rectangles;
    * n_var : (int) ==> number of values of variable var that will be given for each rectangle;
    
    * var : (list) ==> Values related to a property or physical greatness of the layers on the subsurface. These values have to be given in a specific form that will depend of how many rectangles will be painted and how many values of var will be given for each rectangle. Var is a list that contains others list. The number of the the lists inside var provide the number of rectangles that will be painted. The number of elements of these lists provide the number of values related to a property or physical greatness of the layers on the subsurface of the rectangles (This number is the value of the variable n_var that should be given).
    For example:
       In the cases below, five values of var will be given for each rectangle.
       - If you have just one rectangle to paint: var = [ [13,1,10,31,54] ]
       - If you have two rectangles to paint: var = [ [14,25,51,32,61], [2,4,6,7,1] ]
       - If you have three rectangles to paint: var = [ [14,16,71,10,11], [21,22,52,62,29], [32,13,43,23,39] ]
    
    * name : (string) ==> It refers to the name of the property or physical greatness of the layers on the subsurface that is in study. If the rectangles will be painted, you can provide the name of this property or physical greatness as you can see in the examples below. This name will be show on the colorbar.
    For example:
        - name = 'Depth $(m)$'
        - name = 'Density constrast $(g/cm^3)$';
    
    Outputs:

    """
    plt.figure( figsize=(10,10) )
    
    # metade da espessura horizontal do primeiro retangulo em funcao da posição de observacao (1 prisma por observacao):
    # localizacao do primeiro prisma em funcao da posicao de observacao (1 prisma por observacao):
    xmed = ( xo[0] + xo[1] ) / 2.0
    x_prisma = abs( xo[0] - ( xmed ) ) # esse valor será usado durante a etapa de visualização gráfica para a possibilizar a geração dos retangulos
    
    nobs = len(xo)                
    for i in range ( nobs ):
        for j in range( n_var ): 
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
       
    for i in range (nobs):
        if zo[i] < 0:        
            zp = np.array([ zo[i], zo[i], ref + p[i] , ref + p[i], zo[i] ])
            if i == 0:
                zmin = zo[i] 
                zmax = ref + p[i]
            else:
                if ref + p[i]  > zmax:
                    zmax = ref + p[i]
                if zo[i]  < zmin:
                    zmin = zo[i] 
        else:        
            zp = np.array([ zo[i], zo[i], zo[i] + p[i] , zo[i] + p[i], zo[i] ])
            if i == 0:
                zmin = zo[i] 
                zmax = zo[i] + p[i]
            else:
                if zo[i] + p[i]  > zmax:
                    zmax = zo[i] + p[i]
                if zo[i]  < zmin:
                    zmin = zo[i] 
    
        xp = np.array([xo[i] - x_prisma, xo[i] + x_prisma, xo[i] + x_prisma, xo[i] - x_prisma, xo[i] - x_prisma])
    
        path = Path(np.array([xp,zp]).T)
        patch = PathPatch(path, facecolor='none')
    
        plt.gca().add_patch(patch)
        fs = 18 # font size for the label
        plt.ylabel('Depth $(m)$',fontsize=fs)
    
        var_part = np.array( var[i] ) # Isolando apenas os dados do contraste de densidade que serão usados nesse for  
    
        im = plt.imshow(var_part.reshape(n_var,1), cmap, interpolation="bicubic", 
                    vmin=var_min, vmax=var_max,
                    origin='lower',extent=[min(xp), max(xp), min(zp), max(zp)],
                    aspect="auto", clip_path=patch, clip_on=True)

    cbar = plt.colorbar()
    cbar.ax.set_ylabel(name, fontsize=fs)

    plt.xlim( [ xo[0] - 2 * x_prisma , xo[nobs - 1] + 2 * x_prisma ] ) # retirar esse comando, interfere na visualização
    plt.ylim(zmax + 2, zmin - 2)

    plt.grid()
    plt.show()
        
###########################################################################################################################################
def teste_plots_rectangles(xo, zo, p, ref=None, color1='blue', color2='black', color3='green'):        

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
        z_plot.append( [ zo[i], zo[i], zo[i] + p[i] , zo[i] + p[i], zo[i] ] )
        
    # Visualização gráfica:
    plt.figure( figsize=(10,10) )
    plt.plot(xo,zo,'vr')
    
    plt.plot(xo, zo, color1) # delimitação do relevo
   
    z_lim_sed_emb = []
    for i in range(nobs):
        z_lim_sed_emb.append( z_plot[i][3] )
    
    plt.plot(xo, z_lim_sed_emb, color3) # delimitação do limite sedimento - embasamento
    
    for i in range(nobs):
        plt.plot( x_plot[i], z_plot[i], color2) # visualização dos retângulos
    
    plt.grid()
    plt.xlim( [ xo[0] - 2 * x_prisma , xo[nobs - 1] + 2 * x_prisma ] ) # retirar esse comando, interfere na visualização
    plt.gca().invert_yaxis()
    plt.show()
       ###########################################################################################################################################
#def teste_plots_paint_rectangles(xo, zo, p, ref=0, n_var=None, var=None, name=None, cmap='RdBu_r'):        
        
        
###########################################################################################################################################        