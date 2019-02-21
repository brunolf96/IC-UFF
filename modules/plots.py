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
def plots_rectangles(xo, zo, p, color='black' ):
    """
    Plot a group of rectangles with the same width in a juxtaposed framework following a number of observations.

    Inputs:
    * xo : (1D array) ==>  x coordinate of observations;
    * zo : (1D array) ==>  z coordinate of observations;   
    *  p : (1D array) ==>  difference between top and bottom of a specific rectangle;     
    * color : (string) ==> color of the rectangles 
    Outputs:

    """

    # localizacao do primeiro prisma em funcao da posicao de observacao (1 prisma por observacao):
    xmed = ( xo[0] + xo[1] ) / 2.0
    x_prisma = xo[0] - ( xmed )
    

    # Visualização gráfica:
    plt.figure( figsize=(10,10) )
    plt.plot(xo,zo,'vr')
    nobs = len(xo)
    for i in range(nobs):
        plt.plot((xo[i] - x_prisma, xo[i] + x_prisma, xo[i] + x_prisma, xo[i] - x_prisma, xo[i] - x_prisma), 
                 (zo[i], zo[i], zo[i]+p[i] , zo[i]+p[i], zo[i]), color)

    plt.grid()
    #plt.ylim(min(zo), max(zo) + p[i])
    plt.gca().invert_yaxis()
    plt.show()