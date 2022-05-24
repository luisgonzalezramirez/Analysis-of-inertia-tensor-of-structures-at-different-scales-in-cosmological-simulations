
#Código de lectura de datos y realización de gráficas:
#_______________________________________________________________#

#Importar módulos de libreria de python:

import time
import matplotlib.pyplot as plt
import os
from os import listdir
import shutil
import struct
import numpy as np

start = time.time()

def clearConsole():
    command = 'clear'
    if os.name in ('nt', 'dos'):
        command = 'cls'
    os.system(command)

def Vector_Normalizado(v):
    ''' Devuelve el vector normalizado'''
    v = v / np.linalg.norm(v)
    return v

def Ángulo(v1, v2):
    ''' Devuelve el ángulo entre los vectores v1 y v2 en radianes:

            >>> Ángulo((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> Ángulo((1, 0, 0), (1, 0, 0))
            0.0
            >>> Ángulo((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    '''
    v1_u = Vector_Normalizado(v1)
    v2_u = Vector_Normalizado(v2)
    alpha = np.arccos(abs(np.dot(v1_u, v2_u)))
    return alpha
    
def Componentes_Web_no_print(cgx,cgy,cgz,xweb_array,yweb_array,zweb_array,rmweb_array):
    
    '''
    Función:
    ________________________________________________________________________________________
    Realiza el cálculo de las componentes del tensor de inercia I_ij.
    
    Esta función importa el módulo scipy.linalg para la obtención de autovalores.
    
    Variables:
    ________________________________________________________________________________________
    * cgx, coordenada x del centro de masas.
    * cgy, coordenada y del centro de masas.
    * cgz, coordenada z del centro de masas.
    * xweb_array, array de las componentes x de las estructuras web.
    * yweb_array, array de las componentes y de las estructuras web.
    * zweb_array, array de las componentes z de las estructuras web.
    * rmweb_array, array de las masas de las estructuras web.
    '''
    
    import scipy.linalg as la
    
    I_xx_n_array = []
    I_yy_n_array = []
    I_zz_n_array = []
    I_xy_n_array = []
    I_yz_n_array = []
    I_zx_n_array = []
    
    #I^r_ij = sum_n m_n (delta_ij r^2_n - r_in r_jn) / r^2_n
    
    M = np.sum(rmweb_array)
    
    for j in np.linspace(0,len(xweb_array)-1,len(xweb_array)).astype(int):
        r_n2 = (xweb_array[j] - cgx) ** 2 + (yweb_array[j] - cgy) ** 2 + (zweb_array[j] - cgz)** 2
        r_x = abs(xweb_array[j] - cgx)
        r_y = abs(yweb_array[j] - cgy)
        r_z = abs(zweb_array[j] - cgz)
        m_n = rmweb_array[j]
        
        I_xx_n = np.multiply(m_n,(r_n2 - (r_x **2)))/r_n2
        I_xx_n_array.append(I_xx_n)
        
        I_yy_n = np.multiply(m_n,(r_n2 - (r_y **2)))/r_n2
        I_yy_n_array.append(I_yy_n)
        
        I_zz_n = np.multiply(m_n,(r_n2 - (r_z **2)))/r_n2
        I_zz_n_array.append(I_zz_n)
        
        I_xy_n = np.multiply(m_n,(-r_x * r_y))/r_n2
        I_xy_n_array.append(I_xy_n)
        
        I_yz_n = np.multiply(m_n,(-r_y * r_z))/r_n2
        I_yz_n_array.append(I_yz_n)
        
        I_zx_n = np.multiply(m_n,(-r_z * r_x))/r_n2
        I_zx_n_array.append(I_zx_n)
        
    I_xx = np.sum(I_xx_n_array)
    I_yy = np.sum(I_yy_n_array)
    I_zz = np.sum(I_zz_n_array)
    I_xy = np.sum(I_xy_n_array)
    I_yz = np.sum(I_yz_n_array)
    I_zx = np.sum(I_zx_n_array)
    
    I = np.array([[I_xx,I_xy,I_zx],[I_xy,I_yy,I_yz],[I_zx,I_yz,I_zz]])
    
    evals, evecs = la.eig(I)
    evals = evals.real
    
    evals, evecs = zip(*sorted(zip(evals, evecs)))
    
    a = np.sqrt((np.multiply(5,(evals[1] - evals[0] + evals[2])))/np.multiply(2,M))
    b = np.sqrt((np.multiply(5,(evals[2] - evals[1] + evals[0])))/np.multiply(2,M))
    c = np.sqrt((np.multiply(5,(evals[0] - evals[2] + evals[1])))/np.multiply(2,M))
    
    
    return evals,evecs,a,b,c
    
#Input de los directorios para separar los archivos:
    
print('___________________________________________________________________________________________________________')
directory = input('Introduce dirección local de archivos:\n')

print('___________________________________________________________________________________________________________')
target = input('Introduce dirección local del código:\n')

print('___________________________________________________________________________________________________________')
structure = input('Introduce dirección local de los archivos web:\n')

print('___________________________________________________________________________________________________________')
images_path = input('Introduce dirección local donde guardar gráficos:\n')

print('___________________________________________________________________________________________________________')
satellites = input('Introduce dirección local de las listas de satélites:\n')

frames = listdir(directory)

webs = listdir(structure)

sats = listdir(satellites)

#Copia los archivos del directorio de los satélites en la carpeta local del código:
        
for root, directories, sats in os.walk(satellites, topdown=False):
    for name in sats:
        f = os.path.join(root, name)
        shutil.copy2(f, target)

#Copia los archivos del directorio de la simulación en la carpeta local del código:

for root, directories, frames in os.walk(directory, topdown=False):
    for name in frames:
        f = os.path.join(root, name)
        shutil.copy2(f, target)
      
#Copia los archivos del directorio de las estructuras en la carpeta local del código:
        
for root, directories, webs in os.walk(structure, topdown=False):
    for name in webs:
        f = os.path.join(root, name)
        shutil.copy2(f, target)

#Recordatorio: 'C:/Users/lgrjr/Desktop/UCM/Curso 2021-2022/TFG/Simulación/':

#Se define una lista para los redshift:
redshift_array = []

#Lectura de todas las entradas de la simulación:
for frame in frames:
    
    clearConsole()
    
    f = open(frame, "rb")
    
    #Ficheros para exportar
        #Centros de masas de las estructuras web.
        #Autovalores y semiejes.
    
    #parameter(G=6.67e-8,
     #&          Mpc=3.086e24,
     #&          H0=1e7/Mpc,
     #&          yr=3.16e7,
     #&          kb=1.38e-16,
     #&          mum=1.22*1.672649e-24,
     #&          Msun=2e33)


# Cálculo de la densidad crítica y el tamaño de la caja (en Mpc)
# ----------------------------------------------
#     rhoc=3*H0**2*h100**2/(8.*pi*G)
#     box=box100/h100
# Cálculo de la unidad de longitud (longitud de celda=1)
# ----------------------
#     fl=FLOAT(L)-2.*padding
#     lunit=box*Mpc/fl
# Cálculo de las unidades de masa y tiempo
# -------------------------------
#     munit=(rhoc*Mpc)*(Mpc/Msun)*Mpc*box**3/rmtot*omega0
#     tunit=h0t0/H0/h100
# Unidades de velocidad y energía específica
# ---------------------------------
#     vunit=lunit/tunit
#     eunit=vunit**2

# Unidad de densidad másica
# -----------------
#     dunit=lunit*(lunit/Msun)*(lunit/munit)

#     convL=lunit*fl*atime/Mpc  !Conversion a Mpc
#     convM=munit/1e10          !Conversion a 10^10 M_o
#     convV=vunit*fl*atime/100000. !Conversion a kms^-1
#     convT=Kunit*(fl*atime)**2 !Conversion a K

    #ESTRUCTURAS WEB:
        #web2 - filamento
        #web3 - plano

    #Lectura de parámetros:
    parameters1 = f.read(20)
    l1, out, itime, irun, l2 = struct.unpack(">iiiii", parameters1)
    
    parameters2 = f.read(20)
    l1, out, itstop, nobj, l2 = struct.unpack(">iiiii", parameters2)


    parameters3 = f.read(20)
    l1, out, itdump, ngas, l2 = struct.unpack(">iiiii", parameters3)

    parameters4 = f.read(20)
    l1, out, iout, ndark, l2 = struct.unpack(">iiiii", parameters4)

    n_t=ndark+ngas

    parameters5 = f.read(20)
    l1, out, rtime, L, l2 = struct.unpack(">iiffi", parameters5)

    parameters6 = f.read(20)
    l1, out, atime, intl, l2 = struct.unpack(">iifii", parameters6)
    
    #Redshift:
    z = (1-atime)/atime
    
    redshift_array.append(z)
    
    #Lectura del resto de líneas de parámetros (la mayoría vacías):
    p = range(1,95,1)
    for i in p:
        f.read(20)
    
    r = range(0,nobj,1)         
    
    #Se definen los arrays para las partículas:
    m_array = [] #Masa.
    x_array = [] #Posición x.
    y_array = [] #Posición y.
    z_array = [] #Posición z.
    type_array = [] #Tipo de partícula.
    
    #Se leen los datos de las partículas:
    for i in r:
        data = f.read(40)
        l1,rm,rvx,rvy,rvz,vx,vy,vz,typ,l2 = struct.unpack(">ifffffffii", data)
        m_array.append(rm)
        x_array.append(rvx)
        y_array.append(rvy)
        z_array.append(rvz)
        type_array.append(typ)
        #Se pueden agrupar en una matriz de arrays, pero consume mucha memoria.
    
    #Se define una lista para el número de partículas de la estructura y los satélites:
    web_array = []   
    web_arrayconverted = []
    sat_array = []   
    sat_arrayconverted = []

    #Obtención de partículas de las estructuras web2 y web3, correspondientes a un filamento y un plano:
    for web in webs:
        fweb = open(web, "r")
        #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
        num_par = fweb.readlines()
        #Se ha definido una lista nueva para el número de partículas eliminando la posición 0 que es el valor de la longitud:
        del num_par[0]
        #Se crean un array para los números de las partículas de cada estructura:
        converted_num_par_array = np.array(num_par)
        converted_num_par_array = converted_num_par_array.astype(int)
        #Se obtiene un array de dos arrays con las partículas de las dos estructuras:
        web_arrayconverted.append(converted_num_par_array)
        web_array.append(num_par)
        length = len(num_par)
        rweb = range(len(num_par))

    #Obtención de partículas de los satélites persistentes y no persistentes:
    for sat in sats:
        
        fsat = open(sat, "r")
        #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
        num_par_sat = fsat.readlines()
        #Se crean un array para los números de las partículas de cada estructura:
        converted_num_par_sat_array = np.array(num_par_sat)
        converted_num_par_array = converted_num_par_array.astype(int)
        #Se obtiene un array de dos arrays con las partículas de las dos estructuras:
        sat_arrayconverted.append(converted_num_par_sat_array)
        sat_array.append(num_par_sat)
        length_sat = len(num_par_sat)
        rsat = range(len(num_par_sat))
    
    
    #Se definen arrays para las estructuras web.
    rmweb_array = [] #Masas de estructuras web.
    xweb_array = [] #Posición x estructura web.
    yweb_array = [] #Posición y estructura web.
    zweb_array = [] #Posición z estructura web.
    
    
    i = 0
    n = 2
    
    #Se definen arrays para las partículas de satélites.
    rmsat_array = [] #Masas de las partículas de satélites.
    xsat_array = [] #Posición x de las partículas de satélites.
    ysat_array = [] #Posición y de las partículas de satélites.
    zsat_array = [] #Posición z de las partículas de satélites.
    
    
    #Se obtienen los arrays de las estructuras web y satélites a partir de array de dos arrays con los números de partículas:
    for web,sat in zip(web_arrayconverted, sat_arrayconverted):
        for j,k in zip(web,sat):
            mweb = m_array[int(j)]
            rmweb_array.append(mweb)
            xweb = x_array[int(j)]
            xweb_array.append(xweb)
            yweb = y_array[int(j)]
            yweb_array.append(yweb)
            zweb = z_array[int(j)]
            zweb_array.append(zweb)
            msat = m_array[int(k)]
            rmsat_array.append(msat)
            xsat = x_array[int(k)]
            xsat_array.append(xsat)
            ysat = y_array[int(k)]
            ysat_array.append(ysat)
            zsat = z_array[int(k)]
            zsat_array.append(zsat)
        
        #Se calculan los centros de masa iniciales:
        cgx = np.sum(np.multiply(xweb_array,rmweb_array))/np.sum(rmweb_array)
        cgy = np.sum(np.multiply(yweb_array,rmweb_array))/np.sum(rmweb_array)
        cgz = np.sum(np.multiply(zweb_array,rmweb_array))/np.sum(rmweb_array)
        
        #Se calculan los centros de masa:
        cgx_sat = np.sum(np.multiply(xsat_array,rmsat_array))/np.sum(rmsat_array)
        cgy_sat = np.sum(np.multiply(ysat_array,rmsat_array))/np.sum(rmsat_array)
        cgz_sat = np.sum(np.multiply(zsat_array,rmsat_array))/np.sum(rmsat_array)
        
        #Se dibujan las estructuras web:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xweb_array, yweb_array, zweb_array, alpha=1, marker='.', linewidths=0.5, s=rmweb_array);
        ax.scatter(xsat_array, ysat_array, zsat_array, alpha=1, marker='.', color = 'r', linewidths=0.5, s=rmsat_array);
        plt.xlim(0,1)
        plt.ylim(0,1)
        ax.set_zlim3d(0,1)
        #y el centro de masas inicial:
        ax.scatter(cgx, cgy, cgz, color='k', marker='+', s=1e3,linewidths=0.5);
        ax.scatter(cgx_sat, cgy_sat, cgz_sat, color='r', marker='+', s=1e3,linewidths=0.5);
        #junto con los vectores:
        evals,evecs,a,b,c = Componentes_Web_no_print(cgx,cgy,cgz,xweb_array,yweb_array,zweb_array,rmweb_array)
        evals_sat,evecs_sat,a_sat,b_sat,c_sat = Componentes_Web_no_print(cgx_sat,cgy_sat,cgz_sat,xsat_array,ysat_array,zsat_array,rmsat_array)
        ax.quiver(cgx, cgy, cgz, evecs[0]/(3 * np.linalg.norm(evecs[0])), evecs[1]/(3 * np.linalg.norm(evecs[1])), evecs[2]/(3 * np.linalg.norm(evecs[2])), color='k', linewidth=0.5, arrow_length_ratio=0.1,)
        ax.quiver(cgx_sat, cgy_sat, cgz_sat, evecs_sat[0]/(3 * np.linalg.norm(evecs_sat[0])), evecs_sat[1]/(3 * np.linalg.norm(evecs_sat[1])), evecs_sat[2]/(3 * np.linalg.norm(evecs_sat[2])), color='m', linewidth=0.5, arrow_length_ratio=0.1,)
        plt.title('Estructura {} vs {}, z = {:.2f}'.format(webs[i],sats[i], z))
        plt.savefig(images_path + 'Estructura web{:.2f} vs {}, frame {},{}.png'.format(n, sats[i], itime, irun), dpi=300)
        plt.close()

        i = i + 1
        n = n + 1

#Se define una lista para los redshift:
redshift_array = []

#Lectura de todas las entradas de la simulación:
for frame in frames:
    
    clearConsole()
    
    f = open(frame, "rb")
    
    #Ficheros para exportar
        #Centros de masas de las estructuras web.
        #Autovalores y semiejes.
    
    #parameter(G=6.67e-8,
     #&          Mpc=3.086e24,
     #&          H0=1e7/Mpc,
     #&          yr=3.16e7,
     #&          kb=1.38e-16,
     #&          mum=1.22*1.672649e-24,
     #&          Msun=2e33)


# Cálculo de la densidad crítica y el tamaño de la caja (en Mpc)
# ----------------------------------------------
#     rhoc=3*H0**2*h100**2/(8.*pi*G)
#     box=box100/h100
# Cálculo de la unidad de longitud (longitud de celda=1)
# ----------------------
#     fl=FLOAT(L)-2.*padding
#     lunit=box*Mpc/fl
# Cálculo de las unidades de masa y tiempo
# -------------------------------
#     munit=(rhoc*Mpc)*(Mpc/Msun)*Mpc*box**3/rmtot*omega0
#     tunit=h0t0/H0/h100
# Unidades de velocidad y energía específica
# ---------------------------------
#     vunit=lunit/tunit
#     eunit=vunit**2

# Unidad de densidad másica
# -----------------
#     dunit=lunit*(lunit/Msun)*(lunit/munit)

#     convL=lunit*fl*atime/Mpc  !Conversion a Mpc
#     convM=munit/1e10          !Conversion a 10^10 M_o
#     convV=vunit*fl*atime/100000. !Conversion a kms^-1
#     convT=Kunit*(fl*atime)**2 !Conversion a K

    #ESTRUCTURAS WEB:
        #web2 - filamento
        #web3 - plano

    #Lectura de parámetros:
    parameters1 = f.read(20)
    l1, out, itime, irun, l2 = struct.unpack(">iiiii", parameters1)
    
    parameters2 = f.read(20)
    l1, out, itstop, nobj, l2 = struct.unpack(">iiiii", parameters2)


    parameters3 = f.read(20)
    l1, out, itdump, ngas, l2 = struct.unpack(">iiiii", parameters3)

    parameters4 = f.read(20)
    l1, out, iout, ndark, l2 = struct.unpack(">iiiii", parameters4)

    n_t=ndark+ngas

    parameters5 = f.read(20)
    l1, out, rtime, L, l2 = struct.unpack(">iiffi", parameters5)

    parameters6 = f.read(20)
    l1, out, atime, intl, l2 = struct.unpack(">iifii", parameters6)
    
    #Redshift:
    z = (1-atime)/atime
    
    redshift_array.append(z)
    
    #Lectura del resto de líneas de parámetros (la mayoría vacías):
    p = range(1,95,1)
    for i in p:
        f.read(20)
    
    r = range(0,nobj,1)       
    
    #Se definen los arrays para las partículas:
    m_array = [] #Masa.
    x_array = [] #Posición x.
    y_array = [] #Posición y.
    z_array = [] #Posición z.
    type_array = [] #Tipo de partícula.
    
    #Se leen los datos de las partículas:
    for i in r:
        data = f.read(40)
        l1,rm,rvx,rvy,rvz,vx,vy,vz,typ,l2 = struct.unpack(">ifffffffii", data)
        m_array.append(rm)
        x_array.append(rvx)
        y_array.append(rvy)
        z_array.append(rvz)
        type_array.append(typ)
        #Se pueden agrupar en una matriz de arrays, pero consume mucha memoria.
    
    #Se define una lista para el número de partículas de la estructura y los satélites:
    web_array = []   
    web_arrayconverted = []
    sat_array = []   
    sat_arrayconverted = []

    #Obtención de partículas de las estructuras web2 y web3, correspondientes a un filamento y un plano:
    for web in webs:
        fweb = open(web, "r")
        #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
        num_par = fweb.readlines()
        #Se ha definido una lista nueva para el número de partículas eliminando la posición 0 que es el valor de la longitud:
        del num_par[0]
        #Se crean un array para los números de las partículas de cada estructura:
        converted_num_par_array = np.array(num_par)
        converted_num_par_array = converted_num_par_array.astype(int)
        #Se obtiene un array de dos arrays con las partículas de las dos estructuras:
        web_arrayconverted.append(converted_num_par_array)
        web_array.append(num_par)
        length = len(num_par)
        rweb = range(len(num_par))

    #Obtención de partículas de los satélites persistentes y no persistentes:
    for sat in sats:
        
        fsat = open(sat, "r")
        #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
        num_par_sat = fsat.readlines()
        #Se crean un array para los números de las partículas de cada estructura:
        converted_num_par_sat_array = np.array(num_par_sat)
        converted_num_par_array = converted_num_par_array.astype(int)
        #Se obtiene un array de dos arrays con las partículas de las dos estructuras:
        sat_arrayconverted.append(converted_num_par_sat_array)
        sat_array.append(num_par_sat)
        length_sat = len(num_par_sat)
        rsat = range(len(num_par_sat))
    
    
    
    
    
    i = 0
    n = 2
    
    
    element1 = sat_arrayconverted[0]
    element0 = sat_arrayconverted[1]
    
    sat_arrayconverted_order = [element0, element1]
    
    #Se obtienen los arrays de las estructuras web y satélites a partir de array de dos arrays con los números de partículas:
    for web,sat in zip(web_arrayconverted, sat_arrayconverted_order):
        #Se definen arrays para las estructuras web.
        rmweb_array = [] #Masas de estructuras web.
        xweb_array = [] #Posición x estructura web.
        yweb_array = [] #Posición y estructura web.
        zweb_array = [] #Posición z estructura web.
        #Se definen arrays para las partículas de satélites.
        rmsat_array = [] #Masas de las partículas de satélites.
        xsat_array = [] #Posición x de las partículas de satélites.
        ysat_array = [] #Posición y de las partículas de satélites.
        zsat_array = [] #Posición z de las partículas de satélites.
        for j,k in zip(web,sat):
            mweb = m_array[int(j)]
            rmweb_array.append(mweb)
            xweb = x_array[int(j)]
            xweb_array.append(xweb)
            yweb = y_array[int(j)]
            yweb_array.append(yweb)
            zweb = z_array[int(j)]
            zweb_array.append(zweb)
            msat = m_array[int(k)]
            rmsat_array.append(msat)
            xsat = x_array[int(k)]
            xsat_array.append(xsat)
            ysat = y_array[int(k)]
            ysat_array.append(ysat)
            zsat = z_array[int(k)]
            zsat_array.append(zsat)
        
        #Se calculan los centros de masa iniciales:
        cgx = np.sum(np.multiply(xweb_array,rmweb_array))/np.sum(rmweb_array)
        cgy = np.sum(np.multiply(yweb_array,rmweb_array))/np.sum(rmweb_array)
        cgz = np.sum(np.multiply(zweb_array,rmweb_array))/np.sum(rmweb_array)
        
        #Se calculan los centros de masa:
        cgx_sat = np.sum(np.multiply(xsat_array,rmsat_array))/np.sum(rmsat_array)
        cgy_sat = np.sum(np.multiply(ysat_array,rmsat_array))/np.sum(rmsat_array)
        cgz_sat = np.sum(np.multiply(zsat_array,rmsat_array))/np.sum(rmsat_array)
        
        #Se dibujan las estructuras web:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xweb_array, yweb_array, zweb_array, alpha=1, marker='.', linewidths=0.5, s=rmweb_array);
        ax.scatter(xsat_array, ysat_array, zsat_array, alpha=1, marker='.', color = 'r', linewidths=0.5, s=rmsat_array);
        plt.xlim(0,1)
        plt.ylim(0,1)
        ax.set_zlim3d(0,1)
        #y el centro de masas inicial:
        ax.scatter(cgx, cgy, cgz, color='k', marker='+', s=1e3,linewidths=0.5);
        ax.scatter(cgx_sat, cgy_sat, cgz_sat, color='r', marker='+', s=1e3,linewidths=0.5);
        #junto con los vectores:
        evals,evecs,a,b,c = Componentes_Web_no_print(cgx,cgy,cgz,xweb_array,yweb_array,zweb_array,rmweb_array)
        evals_sat,evecs_sat,a_sat,b_sat,c_sat = Componentes_Web_no_print(cgx_sat,cgy_sat,cgz_sat,xsat_array,ysat_array,zsat_array,rmsat_array)
        ax.quiver(cgx, cgy, cgz, evecs[0]/(3 * np.linalg.norm(evecs[0])), evecs[1]/(3 * np.linalg.norm(evecs[1])), evecs[2]/(3 * np.linalg.norm(evecs[2])), color='k', linewidth=0.5, arrow_length_ratio=0.1,)
        ax.quiver(cgx_sat, cgy_sat, cgz_sat, evecs_sat[0]/(3 * np.linalg.norm(evecs_sat[0])), evecs_sat[1]/(3 * np.linalg.norm(evecs_sat[1])), evecs_sat[2]/(3 * np.linalg.norm(evecs_sat[2])), color='m', linewidth=0.5, arrow_length_ratio=0.1,)
        plt.title('Estructura {} vs {}, z = {:.2f}'.format(webs[i],sats[1 - i], z))
        plt.savefig(images_path + 'Estructura web{:.2f} vs {}, frame {},{}.png'.format(n, sats[1 - i], itime, irun), dpi=300)
        plt.close()

        i = i + 1
        n = n + 1
    
end = time.time()
print('___________________________________________________________________________________________________________')
print("Tiempo de ejecución (s):", end-start)
print("En minutos (min):", (end-start)/60)
print('___________________________________________________________________________________________________________')


