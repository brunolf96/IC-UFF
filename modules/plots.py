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
def plots_rectangles(xo, zo, p, color1='blue', color2='black', color3='green'):        
    """
    Plot a group of rectangles with the same width in a juxtaposed framework following a number of observations.

    Inputs:
    
    * xo : (1D array) ==> x coordinate of observations;
    * zo : (1D array) ==> z coordinate of observations;   
    *  p : (1D array) ==> difference between top and bottom of a specific rectangle that represents the difference between the relief and the interface that divides the sediment layer from the basement layer;
    * color1 : (string) ==> color of the curve that represents the relief; 
    * color2 : (string) ==> color of the rectangles; 
    * color3 : (string) ==> color of the curve that represents the interface that divides the sediment layer from the basement layer;
    
    Outputs:
    
    * x_plot : (list) ==> It is a list that contains other lists where each of them provides the x coordinate of the points that form a specific rectangle that has been plotted;
    * z_plot : (list) ==> It is a list that contains other lists where each of them provides the z coordinate of the points that form a specific rectangle that has been plotted;

    """
    # cálculo do valor que deve ser somado ou subtraido aos pontos de observação para se ter as coordenadas x dos retângulos:
    xmed = ( xo[0] + xo[1] ) / 2.0
    x_prisma = abs( xo[0] - ( xmed ) ) # esse valor será usado durante a etapa de visualização gráfica para a possibilizar a geração dos retangulos
    
    # criando listas que armazenarão as coordenadas usadas nos gráficos produzidos: 
    nobs = len(xo)
    x_plot = []
    z_plot = []
    for i in range(nobs): # construindo listas com todas as coordenadas x e z dos pontos que formam todos os retângulos
        x_plot.append( [ xo[i] - x_prisma, xo[i] + x_prisma, xo[i] + x_prisma, xo[i] - x_prisma, xo[i] - x_prisma ] )
        z_plot.append( [ zo[i], zo[i], zo[i] + p[i] , zo[i] + p[i], zo[i] ] )
        
    # Visualização gráfica:
    plt.figure( figsize=(10,10) )
    plt.plot(xo,zo,'vr') # visualização das observações
    
    plt.plot(xo, zo, color1) # delimitação do relevo
   
    z_lim_sed_emb = []
    for i in range(nobs): # obtenção das coordenadas z dos pontos que formaram a curva que representa o limite sedimento - embasamento
        z_lim_sed_emb.append( z_plot[i][3] )
    
    plt.plot(xo, z_lim_sed_emb, color3) # delimitação do limite sedimento - embasamento
    
    for i in range(nobs):
        plt.plot( x_plot[i], z_plot[i], color2) # visualização dos retângulos
    
    plt.legend(['Observações','Relevo','Sedimento - Embasamento'], loc='upper left', bbox_to_anchor=(1,1)) # legendas
    plt.grid()
    plt.xlim( [ xo[0] - 2 * x_prisma , xo[nobs - 1] + 2 * x_prisma ] ) # retirar esse comando, interfere na visualização
    plt.gca().invert_yaxis()
    plt.show()
    
    return x_plot, z_plot
       ###########################################################################################################################################
def plots_paint_rectangles(xo, zo, p, n_var=None, var=None, name=None, cmap='RdBu_r', color1='blue', color2='green'):        
    """
    Paint a group of rectangles with the same width in a juxtaposed framework.

    Inputs:
    
    * xo : (1D array) ==> x coordinate of observations;
    * zo : (1D array) ==> z coordinate of observations;       
    *  p : (1D array) ==> difference between top and bottom of a specific rectangle that represents the difference between the relief and the interface that divides the layer of sediment from the basement layer;
    * cmap : (string) ==> colormap instance or registered colormap name;
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
    
    * color1 : (string) ==> color of the curve that represents the relief;
    * color2 : (string) ==> color of the curve that represents the interface that divides the sediment layer from the basement layer;
    
    Outputs:

    """      
    plt.figure( figsize=(10,10) )
    
    # cálculo do valor que deve ser somado ou subtraido aos pontos de observação para se ter as coordenadas x dos retângulos:
    xmed = ( xo[0] + xo[1] ) / 2.0
    x_prisma = abs( xo[0] - ( xmed ) ) # esse valor será usado durante a etapa de visualização gráfica para a possibilizar a geração dos retangulos
    
    nobs = len(xo)                
    for i in range ( nobs ):
        # salvando os valores máximo e minimo da propriedade ou grandeza física que será representada pela escala de cor do plt.imshow
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
    
    z_lim_sed_emb = []
    for i in range (nobs):
        zp = np.array([ zo[i], zo[i], zo[i] + p[i] , zo[i] + p[i], zo[i] ])
        # salvando os valores máximo e mínimo de profundidade para definir melhor o intervalo do eixo y
        if i == 0:
            zmin = zo[i] 
            zmax = zo[i] + p[i]
        else:
            if zo[i] + p[i]  > zmax:
                zmax = zo[i] + p[i]
            if zo[i]  < zmin:
                zmin = zo[i] 
        
        z_lim_sed_emb.append( zp[3] ) # obtenção das coordenadas z dos pontos que formaram a curva que representa o limite sedimento - embasamento
    
        xp = np.array([xo[i] - x_prisma, xo[i] + x_prisma, xo[i] + x_prisma, xo[i] - x_prisma, xo[i] - x_prisma])
    
        path = Path(np.array([xp,zp]).T)
        patch = PathPatch(path, facecolor='none')
    
        plt.gca().add_patch(patch)
        fs = 18 # font size for the label
        plt.ylabel('Depth $(m)$',fontsize=fs)
    
        var_part = np.array( var[i] ) # isolando apenas os dados do contraste de densidade que serão usados nesse for  
    
        im = plt.imshow(var_part.reshape(n_var,1), cmap, interpolation="bicubic", 
                    vmin=var_min, vmax=var_max,
                    origin='lower',extent=[min(xp), max(xp), min(zp), max(zp)],
                    aspect="auto", clip_path=patch, clip_on=True)
    
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(name, fontsize=fs)

    plt.xlim( [ xo[0] - 2 * x_prisma , xo[nobs - 1] + 2 * x_prisma ] ) # retirar esse comando, interfere na visualização
    plt.ylim(zmax + 2, zmin - 2)
    
    plt.plot(xo,zo,'vr') # visualização das observações
    plt.plot(xo, zo, color1) # delimitação do relevo
    plt.plot(xo, z_lim_sed_emb, color2) # delimitação do limite sedimento - embasamento
    
    plt.legend(['Observações','Relevo','Sedimento - Embasamento'], loc=9, bbox_to_anchor=(0., .95, 1., .095), ncol=3, mode='expand') # legendas
    plt.grid()
    plt.show()
        
###########################################################################################################################################        