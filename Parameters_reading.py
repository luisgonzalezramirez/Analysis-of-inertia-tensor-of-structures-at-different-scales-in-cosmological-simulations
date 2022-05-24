
#Código de lectura de datos y cálculo de parámetros físicos de la simulación:
#_______________________________________________________________#

#Importar módulos de libreria de python:

import time
import matplotlib.pyplot as plt
import os
from os import listdir
import shutil
import struct
import math
import numpy as np

start = time.time()

def clearConsole():
    command = 'clear'
    if os.name in ('nt', 'dos'):
        command = 'cls'
    os.system(command)

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

#Copia los archivos del directorio de los satélites en la carpeta local del código:
        
for root, directories, sats in os.walk(satellites, topdown=False):
    for name in sats:
        f = os.path.join(root, name)
        shutil.copy2(f, target)


#Se define una lista para los redshift:
redshift_array = []    

#Se define una lista para el número de partículas de la estructura:
web_array = []   
web_arrayconverted = []

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

#Se definen arrays:
array_data = [] #Array para la posibilidad de acumular todo.
cg_array = [] #Centro de masa.
dist_generalarray = [] #Array para distancias entre centro de masas y partículas.
atime_array = [] #Array para la fracción de tiempo.
box100_array = [] #Array con los tamaños de la caja en unidades normalizadas.
h100_array = [] #Array para la constante de Hubble.
L_array = [] #Array para las constantes de longitud.
h0t0_array = [] #Array para las constantes de paso de tiempo.
box_size_array = [] #Array de los tamaños de la caja de la simulación.
omega0_array = [] #Array de las densidades críticas.
total_mass_array = [] #Array de las masas totales.
tunit_array = [] #Array de las unidades de tiempo.

Msun = 2e33
pi = math.pi
G = 6.67e-8
Mpc = 3.086e24
H0 = 1e7/Mpc

#Se definen arrays para las estructuras web para cada z.
rmweb_array_array_array = [] #Masas de estructuras web acumuladas para cada z.
xweb_array_array_array = [] #Posición x estructura web acumuladas para cada z.
yweb_array_array_array = [] #Posición y estructura web acumuladas para cada z.
zweb_array_array_array = [] #Posición z estructura web acumuladas para cada z.
typeweb_array_array_array = [] #Tipos de partículas de estructuras web acumuladas para cada z.

mean_gas_array = []
mean_star_array = []
mean_dm_array = []

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
    
    p = range(0,6,1)
    for i in p:
        f.read(20)
    
    parameters13 = f.read(20)
    l1, out, rad, h100, l2 = struct.unpack(">iiffi", parameters13)
    
    parameters14 = f.read(20)
    l1, out, Esum, box100, l2 = struct.unpack(">iiffi", parameters14)
    
    p = range(0,5,1)
    for i in p:
        f.read(20)
    
    parameters20 = f.read(20)
    l1, out, padding, rmnorm, l2 = struct.unpack(">iiffi", parameters20)
    
    f.read(20)
    
    parameters21 = f.read(20)
    l1, out, Qlost, omega0, l2 = struct.unpack(">iiffi", parameters21)
    
    f.read(20)
    
    parameters24 = f.read(20)
    l1, out, delta, h0t0, l2 = struct.unpack(">iiffi", parameters24)
    
    p = range(0,8,1)
    for i in p:
        f.read(20)
    
    parameters33 = f.read(20)
    l1, out, rmtot, mgrh, l2 = struct.unpack(">iiffi", parameters33)
    
    #Redshift:
    z = (1-atime)/atime
    
    #Densidad crítica y tamaño de la caja:
    rhoc = (3*H0**2)*h100**2/(8.*pi*G)
    box = box100/h100
    
    
    print('Redshift, z =',z)
    print('Fracción de tiempo, t =',atime)
    print('box100 =',box100)
    print('box =',box)
    print('h0t0 =',h0t0)
    print('h100 =',h100)
    print('Omega_0 =',omega0)
    print('Masa total, M =',rmtot)
    print('___________________________________________________________________')
    
    #Unidad de tiempo:
    tunit = h0t0/H0/h100
    tunit_array.append(tunit)
    
    #Unidad de longitud:
    fl=L - 2*padding
    lunit = box*Mpc/fl
    convL = lunit*fl*atime/Mpc
    box_conv = convL
    
    print('Densidad crítica =',rhoc)
    print('Caja =',box_conv)
    print('___________________________________________________________________')
    print('___________________________________________________________________')
    
    redshift_array.append(z)
    atime_array.append(atime)
    box100_array.append(box100)
    L_array.append(L)
    box_size_array.append(box_conv)
    h0t0_array.append(h0t0)
    omega0_array.append(omega0)
    total_mass_array.append(rmtot)
    
    
    
    #Lectura del resto de líneas de parámetros (la mayoría vacías):
    p = range(0,67,1)
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
    
    #Se definen arrays para las estructuras web.
    rmweb_array_array = [] #Masas de estructuras web.
    xweb_array_array = [] #Posición x estructura web.
    yweb_array_array = [] #Posición y estructura web.
    zweb_array_array = [] #Posición z estructura web.
    typeweb_array_array = [] #Tipos de partículas de estructuras web.
    
    m = 0
    
    #Se obtienen los arrays de las estructuras web a partir de array de dos arrays con los números de partículas:
    for web in web_arrayconverted:
        #Se definen arrays para las estructuras web.
        rmweb_array = [] #Masas de estructuras web.
        xweb_array = [] #Posición x estructura web.
        yweb_array = [] #Posición y estructura web.
        zweb_array = [] #Posición z estructura web.
        typeweb_array = [] #Tipos de partículas de estructuras web.
        rmweb_array_array.append(webs[m])
        xweb_array_array.append(webs[m])
        yweb_array_array.append(webs[m])
        zweb_array_array.append(webs[m])
        typeweb_array_array.append(webs[m])
        for j in web:
            mweb = m_array[int(j)]
            rmweb_array.append(mweb)
            xweb = x_array[int(j)]
            xweb_array.append(xweb)
            yweb = y_array[int(j)]
            yweb_array.append(yweb)
            zweb = z_array[int(j)]
            zweb_array.append(zweb)
            typeweb = type_array[int(j)]
            typeweb_array.append(typeweb)
            if typeweb == 0:
                print('..............')
                print('web tienen DM')
                print('..............')
        rmweb_array_array.append(rmweb_array)
        xweb_array_array.append(xweb_array)
        yweb_array_array.append(yweb_array)
        zweb_array_array.append(zweb_array)
        typeweb_array_array.append(typeweb_array)
        m = m + 1
        
    #Posiciones de las partículas de gas:
    xgas_array = []
    ygas_array = []
    zgas_array = []
    #Masas de las partículas de gas:
    mgas_array = []
    #Posiciones de las partículas de material estelar:
    xstar_array = []
    ystar_array = []
    zstar_array = []
    #Masas de las partículas de material estelar:
    mstar_array = []
    #Posiciones de las partículas de materia oscura:
    xdm_array = []
    ydm_array = []
    zdm_array = []
    #Masas de las partículas de materia oscura:
    mdm_array = []
    
    #Dentro del tamaño del array de los tipos de partículas, se seleccionan los tipos para separar las partículas del frame:
    for i in np.linspace(0, len(type_array)-1, len(type_array)).astype(int):
        if type_array[i] == 1:
            xgas = x_array[i]
            xgas_array.append(xgas)
            ygas = y_array[i]
            ygas_array.append(ygas)
            zgas = z_array[i]
            zgas_array.append(zgas)
            mgas = m_array[i]
            mgas_array.append(mgas)
        elif type_array[i] == -1:
            xstar = x_array[i]
            xstar_array.append(xstar)
            ystar = y_array[i]
            ystar_array.append(ystar)
            zstar = z_array[i]
            zstar_array.append(zstar)
            mstar = m_array[i]
            mstar_array.append(mstar)
        elif type_array[i] == 0:
            xdm = x_array[i]
            xdm_array.append(xdm)
            ydm = y_array[i]
            ydm_array.append(ydm)
            zdm = z_array[i]
            zdm_array.append(zdm)
            mdm = m_array[i]
            mdm_array.append(mdm)
    
    n = 0
    
    
    #Se grafica una figura de todas las partículas de la simulación para el redshift correspondiente, añadiendo transparencia
    #para las partículas de materia oscura y de gas:
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    sc1 = ax1.scatter(xgas_array, ygas_array, zgas_array, color='blue', alpha=0.5, marker='.', linewidths=0, s=mgas_array);
    sc2 = ax1.scatter(xstar_array, ystar_array, zstar_array, color='red', marker='.', linewidths=0, s=mstar_array);
    sc3 = ax1.scatter(xdm_array, ydm_array, zdm_array, color='black', alpha=0.05, marker='.', linewidths=0, s=mdm_array);
    plt.xlim(0,1)
    plt.ylim(0,1)
    ax1.set_zlim3d(0,1)
    plt.tick_params(axis="x", labelsize=8)
    plt.tick_params(axis="y", labelsize=8)
    plt.tick_params(axis="z", labelsize=8)
    plt.title('Simulación z = {:.2f}'.format(z))
    plt.savefig(images_path + 'Simulación, frame {},{}.png'.format(itime, irun), dpi=500)
    plt.close()
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    sc1 = ax1.scatter(xgas_array, ygas_array, zgas_array, color='blue', alpha=0.5, marker='.', linewidths=0, s=mgas_array);
    sc2 = ax1.scatter(xstar_array, ystar_array, zstar_array, color='red', marker='.', linewidths=0, s=mstar_array);
    sc3 = ax1.scatter(xdm_array, ydm_array, zdm_array, color='black', alpha=0.05, marker='.', linewidths=0, s=mdm_array);
    plt.xlim(0,1)
    plt.ylim(0,1)
    ax1.set_zlim3d(0,1)
    plt.tick_params(axis="x", labelsize=8)
    plt.tick_params(axis="y", labelsize=8)
    plt.tick_params(axis="z", labelsize=8)
    ax1.view_init(90, 90)
    plt.title('Simulación z = {:.2f}'.format(z))
    plt.savefig(images_path + 'Simulación, frame, cenit {},{}.png'.format(itime, irun), dpi=500)
    plt.close()

end = time.time()
print('___________________________________________________________________________________________________________')
print("Tiempo de ejecución (s):", end-start)
print("En minutos (min):", (end-start)/60)
print('___________________________________________________________________________________________________________')


