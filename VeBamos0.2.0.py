#!/usr/bin/env python3

"""
Created on Fri Jun 11 17:29:00 2021

@author: Camilo Prada Latorre c.prada@uniandes.edu.co
"""
import numpy as np
import os , sys
from scipy.integrate import quad
from scipy.misc import derivative
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import time
from datetime import date
from datetime import datetime

class Molecula:

    '''Esta Clase permite la construccion de los vectores del metodo'''

    #----------------------------------------------------------#
    #-----------------------Constructor------------------------#
    #----------------------------------------------------------#


    def __init__(self, fname):


        #----------------------------------------------------------#
        #------------------------Atributos-------------------------#
        #----------------------------------------------------------#
        self.atomtype= []
        self.coordenadas_atomos = []
        self.volumenMolecula = 0
        self.areaSuperficial = 0
        self.d1 = 0
        self.d2 = []
        self.areaTransversal = 0
        self.Polarizabilidad = 0
        self.MomentoDipolar = 0
        self.MomentoCuadripolar = 0
        self.EnergiaDeDispersion = 0
        self.temp = ''
        self.numhomo = 0
        self.numlumo = 0
        self.HOMO = 0
        self.LUMO = 0
        self.gaphl = 0
        self.dureza = 0
        self.suavidad = 0
        self.potencialquimico = 0
        self.electrofilicidad = 0
        self.vector = np.array([])
        self.cuadrupoloIsotropico = 0
        self.polarizabilidadIsotropica = 0
        self.atomtype = []
        self.coordenadas_atomos = []
        self.nombre = fname
        self.propiedades = []

        #----------------------------------------------------------#
        #------------------------Funciones-------------------------#
        #----------------------------------------------------------#

        self.CargarPropiedadesElectronicas(self.nombre)
        self.CargarCoordenadas(self.nombre)
        self.CalcularD1()
        self.CalcularD2()
        self.CalcularVolumen()
        self.CalcularAreaSuperficial()
        self.GapHomoLumo()
        self.CalcularDurezaGlobal()
        self.CalcularSuavidadGlobal()
        self.CalcularPotencialQuimico()
        self.CalcularElectrofilicidadGlobal()
        self.DefinirVector()


    #----------------------------------------------------------#
    #----------------Leer Coordenadas '.xyz'-------------------#
    #----------------------------------------------------------#


    def CargarCoordenadas(self, nombre):


        xyz = open(inputdir+nombre+'.xyz', 'r')

        coordenadas = xyz.readlines()

        xyz.close()

        coordenadas = coordenadas[2:]


        for line in coordenadas:

            atom , xs, ys, zs = line.split()

            self.atomtype.append(atom)

            self.coordenadas_atomos.append(np.array([float(xs), float(ys), float(zs)]))


    #----------------------------------------------------------#
    #------------------Leer Output de Orca.--------------------#
    #----------------------------------------------------------#

    def CargarPropiedadesElectronicas(self, nombre):

        out = open(inputdir+nombre+'.out', 'r')

        output = out.readlines()

        out.close


        for i in range(len(output)):
            if ('  Orbitals  ' in output[i]):

                a,b,c = output[i].split('...')

                self.numhomo = int(b[0:3])
                self.numlumo = int(b[9:11])

        for i in range(len(output)):
            if ('Dispersion correction  ' in output[i]):

                dispenergy = output[i][32:44]


            elif ('Total Dipole Moment    :' in output[i]):

                a, b, c, d, e, f, g = output[i].split()


                self.dipolarx = float(e)

                self.dipolary = float(f)

                self.dipolarz = float(g)

            elif ('Magnitude (a.u.)' in output[i]):

                a,b,c,d = output[i].split()

                self.magnitudDipolar = float(d)


            elif ('QUADRUPOLE MOMENT (A.U.)' in output[i]):

                a, b, c, d = output[i+10].split()

                self.quadrupolex  = float(a)

                self.quadrupoley  = float(b)

                self.quadrupolez  = float(c)

            elif ('Isotropic quadrupole' in output[i]):

                a,b,c,d = output[i].split()

                self.cuadrupoloIsotropico = float(d)


            elif ('THE POLARIZABILITY TENSOR' in output[i]):

                a, b, c = output[i+8].split()

                self.polarizabilidadx = float(a)

                self.polarizabilidady = float(b)

                self.polarizabilidadz = float(c)

            elif ('Isotropic polarizability' in output[i]):

                a,b,c,d = output[i].split()

                self.polarizabilidadIsotropica = float(d)


            elif(' %d   2.' %float(self.numhomo) in output[i]):

                t,y,u,p = output[i].split()

                self.HOMO = float(u)

            elif(' %d   1.' %float(self.numhomo) in output[i]):

                t,y,u,p = output[i].split()

                self.HOMO = float(u)

            elif(' %d   0.' %float(self.numlumo) in output[i]):

                t,y,u,p = output[i].split()

                self.LUMO = float(u)



        self.EnergiaDeDispersion = float(dispenergy)



  #---------------------------------------------------------------#
  #---=========================================================---#
  #---||-----------------------------------------------------||---#
  #---||---Calculo de propiedades topologicas moleculares.---||---#
  #---||-----------------------------------------------------||---#
  #---=========================================================---#
  #---------------------------------------------------------------#


    #----------------------------------------------------------#
    #-------Calcular Distancia entre puntos mas distantes------#
    #----------------------------------------------------------#

    def CalcularD1(self):

        a = 0

        for i in range(len(self.coordenadas_atomos)):
            j = i
            while (j <= len(self.coordenadas_atomos)-1):   #Recorrer solo media mariz

                sumcuadrados = np.sum( (self.coordenadas_atomos[i] - self.coordenadas_atomos[j])**2 )

                distancia = np.sqrt(sumcuadrados)

                if(a < distancia):
                    a = distancia
                    self.indice1 = i
                    self.indice2 = j

                j = j+1

        self.d1 = a



    #----------------------------------------------------------#
    #-------Calcular Distancias perpendiculares al eje D1.-----#
    #----------------------------------------------------------#

    def CalcularD2(self):

        puntoQ = self.coordenadas_atomos[self.indice1]

        puntoR = self.coordenadas_atomos[self.indice2]

        vecQR = puntoQ - puntoR

        vecQP = []

        a = 0

        for i in range(len(self.coordenadas_atomos)):

            vecQP.append( puntoQ - self.coordenadas_atomos[i])

            w = vecQP[i] - (vecQR * np.dot(vecQP[i], vecQR) / np.dot(vecQR, vecQR))


            if (np.linalg.norm(w) != 0):
                self.d2.append(np.linalg.norm(w))



    #----------------------------------------------------------#
    #-------------Calcular Volumen de la Molecula.-------------#
    #----------------------------------------------------------#



    def CalcularVolumen(self):


#         xx = sp.Symbol('xx')

#         self.ecuacionvolumen = a[0]* (xx**10) + a[1]* (xx**9) + a[2]* (xx**8) + a[3]* (xx**7) + a[4]* (xx**6) + a[5]* (xx**5) + a[6]* (xx**4) + a[7]* (xx**3) + a[8]* (xx**2) + a[9]*xx + a[10]


#         self.volumenMolecula = (sp.integrate(self.ecuacionvolumen**2, (xx,0, len(molepruebas.d2) ))) *np.pi #REVISAR ESTOOOO


        def EcuacionVolumen(xx):

            x = np.linspace(0, len(self.d2), len(self.d2))

            a = np.polyfit(x,self.d2,10)

            return a[0]* (xx**10) + a[1]* (xx**9) + a[2]* (xx**8) + a[3]* (xx**7) + a[4]* (xx**6) + a[5]* (xx**5) + a[6]* (xx**4) + a[7]* (xx**3) + a[8]* (xx**2) + a[9]*xx + a[10]


        def EcuacionVolumen2(xx):

            return (EcuacionVolumen(xx))**2


        if(len(self.coordenadas_atomos) > 4):

            respuesta, error = quad(EcuacionVolumen2, 0, len(self.d2))

            self.volumenMolecula = respuesta * np.pi

        else:

            self.volumenMolecula = (self.d1**3)*(np.pi/6)


    #----------------------------------------------------------#
    #----------------Calcular Area Superficial.----------------#
    #----------------------------------------------------------#

    def CalcularAreaSuperficial(self):

        def EcuacionVolumen(xx):

            x = np.linspace(0, len(self.d2), len(self.d2))

            a = np.polyfit(x,self.d2,10)

#         y = a[0]* (x**10) + a[1]* (x**9) + a[2]* (x**8) + a[3]* (x**7) + a[4]* (x**6) + a[5]* (x**5) + a[6]* (x**4) + a[7]* (x**3) + a[8]* (x**2) + a[9]*x + a[10]


#         plt.figure()
#         plt.plot(x,self.d2, label='distancias D2')
#         plt.plot(x,y, label = 'Curva ajustada')
#         plt.ylabel('D2 (angstroms)')
#         plt.xlabel('Atomos (numeracion) ')
#         plt.title('Funcion ajustada para el volumen.')
#         plt.legend()
#         plt.savefig( './diagramas/'+molepruebas.nombre+'_volumen.png')

            return a[0]* (xx**10) + a[1]* (xx**9) + a[2]* (xx**8) + a[3]* (xx**7) + a[4]* (xx**6) + a[5]* (xx**5) + a[6]* (xx**4) + a[7]* (xx**3) + a[8]* (xx**2) + a[9]*xx + a[10]

        def AreaSuperficial(xx):

            return EcuacionVolumen(xx) * np.sqrt(1 + derivative(EcuacionVolumen,xx,dx=1e-4)**2)

        if(len(self.coordenadas_atomos) > 4):

            tapa1 = 2 * np.pi * self.d2[0]

            tapa2 = 2 * np.pi * self.d2[-1]

            respuesta, error = quad(AreaSuperficial, 0, len(self.d2))

            self.areaSuperficial = respuesta + tapa1 + tapa2

        else:

            self.areaSuperficial = 4 * np.pi * (self.d1/2)**2


  #---------------------------------------------------------------#
  #---=========================================================---#
  #---||-----------------------------------------------------||---#
  #---||---Calculo de propiedades electronicas moleculares---||---#
  #---||-----------------------------------------------------||---#
  #---=========================================================---#
  #---------------------------------------------------------------#


    #----------------------------------------------------------#
    #----------------.Calcular Gap HOMO-LUMO.------------------#
    #----------------------------------------------------------#


    def GapHomoLumo(self):


        self.gaphl = abs(self.HOMO - self.LUMO)


    #----------------------------------------------------------#
    #-----------------Calcular Dureza Global.------------------#
    #----------------------------------------------------------#

    def CalcularDurezaGlobal(self):

        self.dureza = (1/2)*(self.LUMO + self.HOMO)


    #----------------------------------------------------------#
    #----------------Calcular Suavidad Global.-----------------#
    #----------------------------------------------------------#

    def CalcularSuavidadGlobal(self):

        self.suavidad = 1/(2*self.dureza)

    #----------------------------------------------------------#
    #---------------Calcular Potencial Quimico.----------------#
    #----------------------------------------------------------#

    def CalcularPotencialQuimico(self):


        self.potencialquimico = (1/2)*(self.LUMO - self.HOMO)


    #----------------------------------------------------------#
    #------------Calcular Electrofilicidad Global.-------------#
    #----------------------------------------------------------#

    def CalcularElectrofilicidadGlobal(self):


        self.electrofilicidad = (self.potencialquimico**2)/(2*self.dureza)


  #---------------------------------------------------------------#
  #---=========================================================---#
  #---||-----------------------------------------------------||---#
  #---||---Vamo a hacer el vectorcito puej mi doooon. jeje---||---#
  #---||-----------------------------------------------------||---#
  #---=========================================================---#
  #---------------------------------------------------------------#

    def DefinirVector(self):
#        self.vector = np.array([self.polarizabilidadx, self.polarizabilidady, self.polarizabilidadz, self.volumenMolecula])

        self.vector = np.array([self.dipolarx, self.dipolary, self.dipolarz, self.magnitudDipolar, self.EnergiaDeDispersion, self.d1, self.volumenMolecula, self.areaSuperficial, self.quadrupolex, self.quadrupoley, self.quadrupolez, self.cuadrupoloIsotropico, self.polarizabilidadx, self.polarizabilidady, self.polarizabilidadz, self.polarizabilidadIsotropica, self.HOMO, self.LUMO, max(self.d2), self.gaphl, self.dureza, self.suavidad, self.potencialquimico, self.electrofilicidad])



#Día actual
today = date.today()

#Fecha actual
now1 = datetime.now()


# Imágen de bienvenida VeBaMoS.

print("""
  ╔╗  ╔╗ ╔══╗   ╔═╗╔═╗  ╔═══╦╗         |-|   *
  ║╚╗╔╝║ ║╔╗║   ║║╚╝║║  ║╔═╗║║         |-|  _    *  __
  ╚╗║║╔╩═╣╚╝╚╦══╣╔╗╔╗╠══╣╚══╣║         |-|  |  *    |/'
   ║╚╝║║═╣╔═╗║╔╗║║║║║║╔╗╠══╗╠╝        /---\ |~*~~~o~|
   ╚╗╔╣║═╣╚═╝║╔╗║║║║║║╚╝║╚═╝╠╗       /-----\|  O o *|
    ╚╝╚══╩═══╩╝╚╩╝╚╝╚╩══╩═══╩╝      /_______|o___O__|
Un código de similaridad molecular.
""")

# Argumentos y flags para el uso del programa.


parser = argparse.ArgumentParser(description='Vector Based Molecular Similarity (VeBaMoS) es un software de similaridad molecular que permite comparar grupos de moléculas utilizando una descripción vectorial de las moléculas con descriptores mecanico-cuánticos que permiten conseguir información química de las moléculas estudiadas.')


# Arreglar lo del input.
# Arreglar que los argumentos no sean opcionales. (los que no sean)

parser.add_argument('-i', '--input', dest='input', metavar=' ', default='./', help='Dirección del directorio que contiene los outputs y los xyz de ORCA para evaluar la similaridad.')
parser.add_argument('-o', '--output', dest='output', metavar=' ', default='./', help='Nombre del archivo de salida de VeBaMoS, se guardará en la misma carpeta indicada en el input.')
parser.add_argument('-v', '--vectores', dest='vectores', action='store_true', help='Imprimir la matriz de vectores en el output.')
parser.add_argument('-g', '--grafica', dest='grafica', action='store_true', help='Imprimir la gráfica de actividad vs similaridad.')
args = parser.parse_args()

a = 0

for i in range(len(args.input)):
    if(args.input[i]== '/'):
        a = i

a = a+1

inputdir= args.input[0:a]

inputfile= args.input[a:]



#-----------------------------------------------------------------------#
#----------------------Crear los objetos moleculas.---------------------#
#-----------------------------------------------------------------------#


archivos = os.listdir(path=inputdir)

moleculas = []

print('Calculando los vectores moleculares.')

#######################################################################
convension = """
===============================================================================

            ***  Convención de los resultados  ***


Para interpretar los resultados de manera correcta, a continuación se imprime
la manera en la que se representan las diferentes moléculas estudiadas con
VeBaMoS:

"""
num = 0
for archivo in archivos:
    if ('.out' in archivo):
        x = Molecula(archivo[:-4])
        convension += archivo[:-4] + ' es representada por: ' + str(num) + '\n'
        num +=1
        moleculas.append(x)

#######################################################################



p = open(inputdir+args.output,'w')

p.write("""

                ------------------------------
                 ╔╗  ╔╗ ╔══╗   ╔═╗╔═╗  ╔═══╦╗
                 ║╚╗╔╝║ ║╔╗║   ║║╚╝║║  ║╔═╗║║
                 ╚╗║║╔╩═╣╚╝╚╦══╣╔╗╔╗╠══╣╚══╣║
                  ║╚╝║║═╣╔═╗║╔╗║║║║║║╔╗╠══╗╠╝
                  ╚╗╔╣║═╣╚═╝║╔╗║║║║║║╚╝║╚═╝╠╗
                   ╚╝╚══╩═══╩╝╚╩╝╚╝╚╩══╩═══╩╝
              Un código de similaridad molecular
                ------------------------------

        ***********************************************
        |                                             |
        |          Universidad de los Andes           |
        |           Departamento de Química           |
        |                                             |
        |  Grupo de estructura electrónica molecular  |
        |                                             |
        |         Jhon E. Zapata Rivera (Ph.D)        |
        |         Camilo Prada Latorre (B.Sc)         |
        ***********************************************



Queremos agradecer al equipo del Department of theory and spectroscopy en el
Max Planck Institute fuer Kohlenforschung por hacer de Orca un paquete libre,
esto permite obtener la información mecanico-cuántica necesaria para poder
poner en práctica este método de similaridad.


-----------------------------------------------
El directorio del que se se leen los inputs es: {}
-----------------------------------------------


""".format(inputdir))


###################################################################
#       imprimir la convención de las moléculas y los números.
###################################################################


p.write(convension + """
""")


inp = open(args.input, 'r')
inplines = inp.readlines()

p.write("""
================================================================================
   ***Input File:***

""")

propiedades = []

for i in inplines:
    p.write("""║{}""".format(i))


    if('Dipolarx' in i and '#' not in i):
        propiedades.append('Dipolarx')

    if('Dipolary' in i and '#' not in i):
        propiedades.append('Dipolary')

    if('Dipolarz' in i and '#' not in i):
        propiedades.append('Dipolarz')

    if('DipolarT' in i and '#' not in i):
        propiedades.append('DipolarT')

    if('E_Disper' in i and '#' not in i):
        propiedades.append('E_Disper')

    if('Distanc1' in i and '#' not in i):
        propiedades.append('Distanc1')

    if('VolMolec' in i and '#' not in i):
        propiedades.append('VolMolec')

    if('A_Superf' in i and '#' not in i):
        propiedades.append('A_Superf')

    if('Qdrupolx' in i and '#' not in i):
        propiedades.append('Qdrupolx')

    if('Qdrupoly' in i and '#' not in i):
        propiedades.append('Qdrupoly')

    if('Qdrupolz' in i and '#' not in i):
        propiedades.append('Qdrupolz')

    if('QdrupolI' in i and '#' not in i):
        propiedades.append('QdrupolI')

    if('Polarizx' in i and '#' not in i):
        propiedades.append('Polarizx')

    if('Polarizy' in i and '#' not in i):
        propiedades.append('Polarizy')

    if('Polarizz' in i and '#' not in i):
        propiedades.append('Polarizz')

    if('PolarizI' in i and '#' not in i):
        propiedades.append('PolarizI')

    if('EnerHOMO' in i and '#' not in i):
        propiedades.append('EnerHOMO')

    if('EnerLUMO' in i and '#' not in i):
        propiedades.append('EnerLUMO')

    if('Distanc2' in i and '#' not in i):
        propiedades.append('Distanc2')

    if('Gap_HOLU' in i and '#' not in i):
        propiedades.append('Gap_HOLU')

    if('Dureza__' in i and '#' not in i):
        propiedades.append('Dureza__')

    if('Suavidad' in i and '#' not in i):
        propiedades.append('Suavidad')

    if('PQuimico' in i and '#' not in i):
        propiedades.append('PQuimico')

    if('Elfilici' in i and '#' not in i):
        propiedades.append('Elfilici')

toolbar_width = 40

# setup toolbar
sys.stdout.write("[%s]" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

for m in moleculas:
    m.propiedades = propiedades
    sys.stdout.write("-")
    sys.stdout.flush()

sys.stdout.write("]\n") # this ends the progress bar


#-----------------------------------------------------------------------#
#------------------------- Crear los Vectores --------------------------#
#-----------------------------------------------------------------------#

vectores = []

for i in range(len(moleculas)):
    vectores.append(moleculas[i].vector)


#-----------------------------------------------------------------------#
#---------------------- Normalizar los vectores. -----------------------#
#-----------------------------------------------------------------------#


print('Normalizando los vectores.')

vecnorm = []


#######################################################################
toolbar_width = 40

# setup toolbar
sys.stdout.write("[%s]" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

for j in range(len(vectores[i])):
    m=0
    for i in range(len(vectores)):


        m = m + vectores[i][j]




    vecnorm.append(m)

#######################################################################
    sys.stdout.write("-")
    sys.stdout.flush()

sys.stdout.write("]\n") # this ends the progress bar



vecnorm = np.array(vecnorm)




for i in range(len(vectores)):

        vectores[i] = vectores[i]/vecnorm

#-----------------------------------------------------------------------#
#------------------------- SimilVec en Numpy. --------------------------#
#-----------------------------------------------------------------------#

similaridades= []

print('Calculando la similaridad molecular.')



#######################################################################
toolbar_width = 40

# setup toolbar
sys.stdout.write("[%s]" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

for i in range(len(moleculas)):
    for j in range(len(moleculas)):

        alpha = abs(np.linalg.norm(vectores[i]) - np.linalg.norm(vectores[j]))

        v1 = vectores[i]/np.linalg.norm(vectores[i])

        v2 = vectores[j]/np.linalg.norm(vectores[j])

#         print(np.linalg.norm(vectores[i]), i)

        cosbeta = np.dot(v1, v2)

        beta = np.arccos(np.clip(cosbeta, -1, 1))

        z = alpha + beta

        s = 1/np.exp(z*np.log10(len(moleculas)))


        similaridades.append(s)

#######################################################################

    sys.stdout.write("-")
    sys.stdout.flush()

sys.stdout.write("]\n") # this ends the progress bar



similaridades = np.array(similaridades)

similaridades = similaridades.reshape(len(moleculas),len(moleculas))



p.write("""


================================================================================


                      ######################################
                      #   Matriz de Similaridad Molecular  #
                      ######################################

""".format(inputdir))

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#######################################
#   Hay que corregir los espacios antes del número de fila para los demás 'if'  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#######################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i = 0

while (i<= len(similaridades)-5):
    p.write("""
 {: >14} {: >13} {: >13} {: >13} {: >13}
          --------      --------      --------      --------      --------""".format(i, i+1, i+2, i+3, i+4))

    for j in range(len(similaridades)):

        if(len(str(j)) == 1):
            p.write("""
   {} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j], similaridades[i+2][j], similaridades[i+3][j], similaridades[i+4][j]))

        elif(len(str(j)) == 2):   # Para que todo quede en columna, toca poner un espacio menos antes del primer {}, para tenerlo en cuenta si toca programarlo para centenares y millares de moléculas.

            p.write("""
  {} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j], similaridades[i+2][j], similaridades[i+3][j], similaridades[i+4][j]))

        elif(len(str(j)) == 3):

            p.write("""
 {} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j], similaridades[i+2][j], similaridades[i+3][j], similaridades[i+4][j]))

        elif(len(str(j)) == 4):

            p.write("""
{} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j], similaridades[i+2][j], similaridades[i+3][j], similaridades[i+4][j]))

    p.write("""\n""")
    i = i+5

res = (len(similaridades)/5 - int(len(similaridades)/5))*5

res = int(round(res,1))

if(res == 1):

    p.write("""
 {: >14}
          --------""".format(i))

    for j in range(len(similaridades)):

        p.write("""
   {} {: >13.6f}""".format(j, similaridades[i][j]))

    p.write("""\n""")

if(res == 2):
    p.write("""
 {: >14} {: >13}
          --------      --------""".format(i, i+1))

    for j in range(len(similaridades)):

        p.write("""
   {} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j]))

    p.write("""\n""")


if(res == 3):
    p.write("""
 {: >14} {: >13} {: >13}
          --------      --------      --------""".format(i, i+1, i+2))

    for j in range(len(similaridades)):

        p.write("""
   {} {: >13.6f} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j], similaridades[i+2][j]))

    p.write("""\n""")

if(res == 4):

    p.write("""
 {: >14} {: >13} {: >13} {: >13}
          --------      --------      --------      --------""".format(i, i+1, i+2, i+3))

    for j in range(len(similaridades)):

        if(len(str(j)) == 1):
            p.write("""
   {} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j], similaridades[i+2][j], similaridades[i+3][j]))

        if(len(str(j)) == 2):
            p.write("""
  {} {: >13.6f} {: >13.6f} {: >13.6f} {: >13.6f}""".format(j, similaridades[i][j], similaridades[i+1][j], similaridades[i+2][j], similaridades[i+3][j]))

    p.write("""\n""")


############################################################################################
#
#           Impresión de los vectores, con el módulo activado.
#
############################################################################################


if(args.vectores):

    p.write("""


--------------------------------------------------
Se activó la impresión de los vectores moleculares
--------------------------------------------------


================================================================================


                      #####################################
                      #        Vectores Moleculares       #
                      #####################################

    """)

    # propiedades = ['Dipolarx', 'Dipolary', 'Dipolarz', 'DipolarT', 'E_Disper',
    #                 'Distanc1', 'VolMolec', 'A_Superf', 'Qdrupolx', 'Qdrupoly',
    #                 'Qdrupolz', 'QdrupolI', 'Polarizx', 'Polarizy', 'Polarizz',
    #                 'PolarizI', 'EnerHOMO', 'EnerLUMO', 'Distanc2', 'Gap_HOLU',
    #                 'Dureza__', 'Suavidad', 'PQuimico', 'Elfilici']

    i = 0

    while (i<= len(moleculas)-5):
        p.write("""
    {: >14} {: >12} {: >12} {: >12} {: >12}
              --------     --------     --------     --------     --------""".format(i, i+1, i+2, i+3, i+4))

        for j in range(len(propiedades)):

            p.write("""
{} {: >13.6f} {: >12.6f} {: >12.6f} {: >12.6f} {: >12.6f}""".format(propiedades[j], vectores[i][j], vectores[i+1][j], vectores[i+2][j], vectores[i+3][j], vectores[i+4][j]))


        p.write('\n')
        i = i+5


    if(res == 1):

        p.write("""
    {: >14}
              --------""".format(i))

        for j in range(len(propiedades)):

            p.write("""
{} {: >13.6f}""".format(propiedades[j], vectores[i][j]))

        p.write("""\n""")

    if(res == 2):
        p.write("""
    {: >14} {: >12}
              --------     --------""".format(i, i+1))

        for j in range(len(propiedades)):

            p.write("""
{} {: >13.6f} {: >12.6f}""".format(propiedades[j], vectores[i][j], vectores[i+1][j]))

        p.write("""\n""")


    if(res == 3):
        p.write("""
    {: >14} {: >12} {: >12}
              --------     --------     --------""".format(i, i+1, i+2))

        for j in range(len(propiedades)):

            p.write("""
{} {: >13.6f} {: >12.6f} {: >12.6f}""".format(propiedades[j], vectores[i][j], vectores[i+1][j], vectores[i+2][j]))

        p.write("""\n""")

    if(res == 4):

        p.write("""
    {: >14} {: >12} {: >12} {: >12}
              --------     --------     --------     --------""".format(i, i+1, i+2, i+3))

        for j in range(len(propiedades)):

            p.write("""
{} {: >13.6f} {: >12.6f} {: >12.6f} {: >12.6f}""".format(propiedades[j], vectores[i][j], vectores[i+1][j], vectores[i+2][j], vectores[i+3][j]))



        p.write("""\n""")




########### HACE FALTA ARREGLAR ESTO PARA QUE EN VERDAD GRAFIQUE LO QUE DEBE GRAFICAR.

if(args.grafica):

    x = np.linspace(0,100,10000)
    sigma = 1/np.exp(abs(x/10))
    sigman = -(1/np.exp(abs(x/10)))

    plt.figure()
    plt.scatter(x,sigma)
    plt.scatter(x,sigman)
    plt.savefig(inputdir+args.output[:-4]+'.png')


    p.write("""


--------------------------------------------------
Se activó la gráfica actividad vs similaridad
--------------------------------------------------


===============================================================================


#####################################
Se guardó la imagen como: {}
#####################################

    """.format(args.output[:-4]+'.png'))



now2 = datetime.now()

p.write("""

================================================================================


Empezó en:  {}

Terminó en: {}

Duró: {}


                **************************************
                *   VeBaMoS terminó con normalidad   *
                **************************************


""".format(now1,now2,str(now2-now1)[0:7]))


p.close()

# Lograr que solo se calculen las propiedades indicadas.


print('Listo.')
