# --------------------------------------------------------------------------------------------------
# Title: Grav-Mag Codes
# Author: Bruno Lima 
# Colaborador: Rodrigo Bijani 
# Description: Source codes for plotting images.
# --------------------------------------------------------------------------------------------------

# Import Python libraries:
import numpy as np
import matplotlib.pyplot as plt
import pylab as py

######################################################################################################
def plots_rectangles(xo, zo, p, cden, color='black' ):
    """
    Plot a group of rectangles with the same width in a juxtaposed framework following a number of observations.

    Inputs:
    * xo : (1D array) ==>  x coordinate of observations;
    * zo : (1D array) ==>  z coordinate of observations;   
    *  p : (1D array) ==>  difference between top and bottom of a specific rectangle;     
    * cden : (1D array) ==> value of density contrast of a specific rectangle;
    * color : (string) ==> color of the rectangles 
    Outputs:

    """

    # metade da espessura horizontal do primeiro retangulo em funcao da posição de observacao (1 prisma por observacao):
    # localizacao do primeiro prisma em funcao da posicao de observacao (1 prisma por observacao):
    xmed = ( xo[0] + xo[1] ) / 2.0
    x_prisma = xo[0] - ( xmed ) # esse valor será usado durante a etapa de visualização gráfica para a possibilizar a geração dos retangulos
    
    # criando listas que armazenarão as coordenadas usadas nos gráficos produzidos: 
    nobs = len(xo)
    x_plot = []
    y_plot = []
    x_sca = []
    y_sca = []
    z_sca = []
    for i in range(nobs):
        x_plot.append( [ xo[i] - x_prisma, xo[i] + x_prisma, xo[i] + x_prisma, xo[i] - x_prisma, xo[i] - x_prisma ] )
        # O valor 0 foi usado como referência
        # valores acima ou abaixo de 0 para z0, simbolizando um excesso ou falta de relevo, respectivamente, em relação a referência
        if zo[i] < 0:
            y_plot.append( [ zo[i], zo[i], 0 + p[i] , 0 + p[i], zo[i] ] )
        else:
            y_plot.append( [ zo[i], zo[i], zo[i] + p[i] , zo[i] + p[i], zo[i] ] )
        
        x_sca.append( xo[i] )
        y_sca.append( ( y_plot[i][0] + y_plot[i][2] ) / 2 ) # 'y_plot[0]' e 'y_plot[2]' são um dos meios que podem ser usados para representar topo e base
        z_sca.append( cden[i] )
    
    # Visualização gráfica:
    plt.figure( figsize=(10,10) )
    plt.plot(xo,zo,'vr')
    
    for i in range(nobs):
        plt.plot( x_plot[i], y_plot[i], color) # visuação dos retângulos
    
    plt.scatter( x_sca, y_sca, c=z_sca, cmap='nipy_spectral', linewidth=10) # visualização dos valores de constrate de densidade
    cbar = plt.colorbar(fraction=0.03, aspect=30, orientation='vertical')
    cbar.set_label('Density Contrast ($g/cm³$)', fontsize=15)
    plt.grid()
    #plt.ylim(min(zo), max(zo) + p[i])
    plt.gca().invert_yaxis()
    plt.show()