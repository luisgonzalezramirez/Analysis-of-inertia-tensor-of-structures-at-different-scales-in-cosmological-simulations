
#Código de lectura de datos y realización de gráficas para centros de masa:
#_______________________________________________________________#

#Importar módulos de libreria de python:

import time
import matplotlib.pyplot as plt
import os
from os import listdir
import shutil
import struct
import numpy as np
import statistics as st

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

#Recordatorio: 'C:/Users/lgrjr/Desktop/UCM/Curso 2021-2022/TFG/Simulación/':

#Se define una lista para los redshift:
redshift_array = []
#Se define una lista para los tiempos:
atime_array = []

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

colors = ['r','b','g','y']

cgx_web = []
cgy_web = []
cgz_web = []

w_n = 0

#Se definen arrays:
array_data = [] #Array para la posibilidad de acumular todo.
cg_array = [] #Centro de masa.
dist_meanarray = [] #Array para distancias medias entre centro de masas y partículas.

cgx_array = []
cgy_array = []
cgz_array = []

cgx_array_con = []
cgy_array_con = []
cgz_array_con = []
cg_array_con = []

n_tot_web_array_array = []

for frame in frames:
    
    clearConsole()
    
    f = open(frame, "rb")
    
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
    
    atime_array.append(atime)
    
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
    
    #Se definen arrays para las estructuras web.
    
    
    n = 0
    
    dist_mean_array = []
    dist_meanarray = []
    
    n_tot_web_array = []
    
    #Se obtienen los arrays de las estructuras web a partir de array de dos arrays con los números de partículas:
    for web in web_arrayconverted:
        
        rmweb_array = [] #Masas de estructuras web.
        xweb_array = [] #Posición x estructura web.
        yweb_array = [] #Posición y estructura web.
        zweb_array = [] #Posición z estructura web.
        
        for j in web:
            mweb = m_array[int(j)]
            rmweb_array.append(mweb)
            xweb = x_array[int(j)]
            xweb_array.append(xweb)
            yweb = y_array[int(j)]
            yweb_array.append(yweb)
            zweb = z_array[int(j)]
            zweb_array.append(zweb)
        
        n_tot_web = len(rmweb_array)
        n_tot_web_array.append(n_tot_web)
        
        #Se calculan los centros de masa iniciales:
        cgx_i = np.sum(np.multiply(xweb_array,rmweb_array))/np.sum(rmweb_array)
        cgy_i = np.sum(np.multiply(yweb_array,rmweb_array))/np.sum(rmweb_array)
        cgz_i = np.sum(np.multiply(zweb_array,rmweb_array))/np.sum(rmweb_array)
        cg_vector = [webs[n],z,cgx_i,cgy_i,cgz_i]
        cg_array.append(cg_vector)
        n = n + 1
        
    dist_generalarray = []
    dist_max_array = []
    dist_generalmaxarray = []
    n_tot_web_array_array.append(n_tot_web_array)
    
    n = 0
    
    #Se calcula las distancias de las partículas de cada estructura al centro de masas inicial: 
    for web in web_arrayconverted:
        dist_mean_array = [webs[n],z]
        dist_max_array = [webs[n],z]
        dist_array = []
        for j in web:
            xweb = x_array[int(j)]
            yweb = y_array[int(j)]
            zweb = z_array[int(j)]
            q1 = np.array([xweb, yweb, zweb])
            q2 = np.array([cgx_i, cgy_i, cgz_i])
            squared_dist = np.sum((q1-q2)**2, axis=0)
            dist = np.sqrt(squared_dist)
            dist_array.append(dist)
        #Se obtiene la distancia media de las partículas de cada estructura web al centro de masas:
        dist_mean = 2*st.mean(dist_array)
        #Se obtiene la distancia máxima de las partículas de cada estructura web al centro de masas:
        dist_max = max(dist_array)
        #Se añade a un vector donde se muestra la estructura web y el redshift correspondiente:
        dist_mean_array.append(dist_mean)
        #Se añade a un vector donde se muestra la estructura web y el redshift correspondiente:
        dist_max_array.append(dist_max)
        #Se crea una matriz con los arrays anteriores para cada redshift y estructura web:
        dist_generalarray.append(dist_mean_array)
        #Se crea una matriz con los arrays anteriores para cada redshift y estructura web:
        dist_generalmaxarray.append(dist_max_array)
        n = n + 1
        
num_i = 0

n_web_array = []

mweb_tot_array_array = []
mweb_tot_corr_array_array = []

#Obtención de partículas de las estructuras web2 y web3, correspondientes a un filamento y un plano:
for web in webs:
    print(webs[w_n])
    w_n = w_n + 1
    
    fweb = open(web, "r")
    #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
    num_par = fweb.readlines()
    #Se ha definido una lista nueva para el número de partículas eliminando la posición 0 que es el valor de la longitud:
    del num_par[0]
    length = len(num_par)
    rweb = range(len(num_par))

    #Se definen arrays:
    array_data = [] #Array para la posibilidad de acumular todo.
    cg_array = [] #Centro de masa.

    cgx_array = []
    cgy_array = []
    cgz_array = []
    
    mweb_tot_corr_array = []
    
    mweb_tot_array = []
    
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
        
        #Se definen arrays para las estructuras web.
        rmweb_array = [] #Masas de estructuras web.
        xweb_array = [] #Posición x estructura web.
        yweb_array = [] #Posición y estructura web.
        zweb_array = [] #Posición z estructura web.
        
        for j in num_par:
            mweb_con = m_array[int(j)]
            rmweb_array.append(mweb_con)
            xweb_con = x_array[int(j)]
            xweb_array.append(xweb_con)
            yweb_con = y_array[int(j)]
            yweb_array.append(yweb_con)
            zweb_con = z_array[int(j)]
            zweb_array.append(zweb_con)
        cgx = np.sum(np.multiply(xweb_array,rmweb_array))/np.sum(rmweb_array)
        cgy = np.sum(np.multiply(yweb_array,rmweb_array))/np.sum(rmweb_array)
        cgz = np.sum(np.multiply(zweb_array,rmweb_array))/np.sum(rmweb_array)
        xweb_array_new = xweb_array
        yweb_array_new = yweb_array
        zweb_array_new = zweb_array
        
        mweb_T = np.sum(rmweb_array)
        
        mweb_tot_array.append(mweb_T)
        
        k = 0
        for (x,y,z) in zip(xweb_array,yweb_array,zweb_array):
            p1 = np.array([x, y, z])
            p2 = np.array([cgx, cgy, cgz])
            squared_dist = np.sum((p1-p2)**2, axis=0)
            dist = np.sqrt(squared_dist)
            #Se determina el límite para que se sitúen las partículas en el 95% del radio máximo:
            if dist >= dist_generalarray[num_i][2]/(dist_generalmaxarray[num_i][2]*0.95):
                del xweb_array_new[k]
                del yweb_array_new[k]
                del zweb_array_new[k]
                del rmweb_array[k]
            k = k + 1
        cgx_con = np.sum(np.multiply(xweb_array_new,rmweb_array))/np.sum(rmweb_array)
        cgy_con = np.sum(np.multiply(yweb_array_new,rmweb_array))/np.sum(rmweb_array)
        cgz_con = np.sum(np.multiply(zweb_array_new,rmweb_array))/np.sum(rmweb_array)
        
        mweb_T_corr = np.sum(rmweb_array)
        
        mweb_tot_corr_array.append(mweb_T_corr)
        
        n_web_array.append(len(xweb_array_new))
        
        cg_vector_con = [z,cgx_con,cgy_con,cgz_con]
        cg_array_con.append(cg_vector_con)
        
        cgx_array.append(cgx_con)
        cgy_array.append(cgy_con)
        cgz_array.append(cgz_con)
        
        #Se dibujan las estructuras web:
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(xweb_array_new, yweb_array_new, zweb_array_new, alpha=1, marker='.', linewidths=0.5, s=rmweb_array);
        #plt.xlim(0,1)
        #plt.ylim(0,1)
        #ax.set_zlim3d(0,1)
        #y el centro de masas inicial:
        #ax.scatter(cgx_con, cgy_con, cgz_con, color='k', marker='+', s=1e3,linewidths=0.5);
        #plt.title('Estructura, conv {}, z = {}'.format(webs[num_i], z))
        #plt.savefig(images_path + 'Estructura,_conv {}, frame {},{}.png'.format(webs[num_i], itime, irun), dpi=300)
        #plt.close()
        
        n = n + 1
        
        #Se obtienen los arrays de las estructuras web a partir de array de dos arrays con los números de partículas:
        #for j in num_par:
        #    mweb = m_array[int(j)]
        #    rmweb_array.append(mweb)
        #    xweb = x_array[int(j)]
        #    xweb_array.append(xweb)
        #    yweb = y_array[int(j)]
        #    yweb_array.append(yweb)
        #    zweb = z_array[int(j)]
        #    zweb_array.append(zweb)
        
        #Se calculan los centros de masa iniciales:
        #cgx = np.sum(np.multiply(xweb_array,rmweb_array))/np.sum(rmweb_array)
        #cgy = np.sum(np.multiply(yweb_array,rmweb_array))/np.sum(rmweb_array)
        #cgz = np.sum(np.multiply(zweb_array,rmweb_array))/np.sum(rmweb_array)
        
        #cgx_array.append(cgx)
        #cgy_array.append(cgy)
        #cgz_array.append(cgz)
        
    num_i = num_i + 1    
    
    mweb_tot_corr_array_array.append(mweb_tot_corr_array)
    
    cgx_web.append(cgx_array)
    cgy_web.append(cgy_array)
    cgz_web.append(cgz_array)

cgx_sat_all = []
cgy_sat_all = []
cgz_sat_all = []

s = 0

#Obtención de partículas de los satélites persistentes y no persistentes:
for sat in sats:
    print(sats[s])
    s = s + 1
    
    fsat = open(sat, "r")
    #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
    num_par_sat = fsat.readlines()
    length_sat = len(num_par_sat)
    rsat = range(len(num_par_sat))

    #Se definen arrays:
    array_data_sat = [] #Array para la posibilidad de acumular todo.
    cg_array_sat = [] #Centro de masa.
    
    cgx_array_sat = []
    cgy_array_sat = []
    cgz_array_sat = []
    
    sat_arrayconverted = []
    sat_array = []
    
    #Obtención de partículas de los satélites, persistentes y no persistentes:
    fsat = open(sat, "r")
    #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
    num_par_sat = fsat.readlines()
    #Se crean un array para los números de las partículas de cada estructura:
    converted_num_par_sat_array = np.array(num_par_sat)
    converted_num_par_sat_array = converted_num_par_sat_array.astype(int)
    #Se obtiene un array de dos arrays con las partículas de las dos grupos de satélites:
    sat_arrayconverted.append(converted_num_par_sat_array)
    sat_array.append(num_par_sat)
    
    #Lectura de todas las entradas de la simulación:
    for frame in frames:
        
        clearConsole()
        
        f = open(frame, "rb")
        
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
        
        #Se definen arrays para las partículas de satélites.
        rmsat_array = [] #Masas de las partículas de satélites.
        xsat_array = [] #Posición x de las partículas de satélites.
        ysat_array = [] #Posición y de las partículas de satélites.
        zsat_array = [] #Posición z de las partículas de satélites.
        
        #Se obtienen los arrays de las partículas de satélites a partir de array de dos arrays con los números de partículas:
        for j in num_par_sat:
            msat = m_array[int(j)]
            rmsat_array.append(msat)
            xsat = x_array[int(j)]
            xsat_array.append(xsat)
            ysat = y_array[int(j)]
            ysat_array.append(ysat)
            zsat = z_array[int(j)]
            zsat_array.append(zsat)
        
        #Se calculan los centros de masa:
        cgx_sat = np.sum(np.multiply(xsat_array,rmsat_array))/np.sum(rmsat_array)
        cgy_sat = np.sum(np.multiply(ysat_array,rmsat_array))/np.sum(rmsat_array)
        cgz_sat = np.sum(np.multiply(zsat_array,rmsat_array))/np.sum(rmsat_array)
        
        cgx_array_sat.append(cgx_sat)
        cgy_array_sat.append(cgy_sat)
        cgz_array_sat.append(cgz_sat)

    cgx_sat_all.append(cgx_array_sat)
    cgy_sat_all.append(cgy_array_sat)
    cgz_sat_all.append(cgz_array_sat)



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#sc = ax.scatter(cgx_web[0], cgy_web[0], cgz_web[0], color=colors[0], alpha=1, marker='.', linewidths=0, s=1);
ax.plot3D(cgx_web[0], cgy_web[0], cgz_web[0], color = colors[0], linestyle = ":", linewidth=0.9, label = "{}".format(webs[0]))
#sc = ax.scatter(cgx_web[1], cgy_web[1], cgz_web[1], color=colors[1], alpha=1, marker='.', linewidths=0, s=1);
ax.plot3D(cgx_web[1], cgy_web[1], cgz_web[1], color = colors[1], linestyle = ":", linewidth=0.9, label = "{}".format(webs[1]))
#sc = ax.scatter(cgx_sat_all[0], cgy_sat_all[0], cgz_sat_all[0], color=colors[2], alpha=1, marker='+', linewidths=0, s=1);
ax.plot3D(cgx_sat_all[0], cgy_sat_all[0], cgz_sat_all[0], color = colors[2], linestyle = ":", linewidth=0.9, label = "{}".format(sats[0]))
#sc = ax.scatter(cgx_sat_all[1], cgy_sat_all[1], cgz_sat_all[1], color=colors[3], alpha=1, marker='+', linewidths=0, s=1);
ax.plot3D(cgx_sat_all[1], cgy_sat_all[1], cgz_sat_all[1], color = colors[3], linestyle = ":", linewidth=0.9, label = "{}".format(sats[1]))
for i in np.linspace(0, len(redshift_array)-1, 6).astype(int):
    ax.text(cgx_sat_all[1][i]-0.015,cgy_sat_all[1][i]-0.015,cgz_sat_all[1][i],'z={}'.format(round(redshift_array[i],2)), size=4, zorder=1,color='k')
    ax.text(cgx_sat_all[0][i]+0.003,cgy_sat_all[0][i]+0.003,cgz_sat_all[0][i],'z={}'.format(round(redshift_array[i],2)), size=4, zorder=1,color='k')         
plt.xlim(0.225,0.35)
plt.ylim(0.425,0.55)
ax.set_zlim3d(0.175,0.3)
ax.set_xlabel(r'$x$', fontsize=10)
ax.set_ylabel(r'$y$', fontsize=10)
ax.set_zlabel(r'$z$', fontsize=10)
ax.tick_params(axis="x", labelsize=8)
ax.tick_params(axis="y", labelsize=8)
ax.tick_params(axis="z", labelsize=8)
plt.legend(loc = "upper left",fontsize=4)
plt.title('Simulación para centros de masa')
plt.savefig(images_path + 'Simulación para centros de masa.png', dpi=500)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#sc = ax.scatter(cgx_web[0], cgy_web[0], cgz_web[0], color=colors[0], alpha=1, marker='.', linewidths=0, s=1);
ax.plot3D(cgx_web[0], cgy_web[0], cgz_web[0], color = colors[0], linestyle = ":", linewidth=0.9, label = "{}".format(webs[0]))
#sc = ax.scatter(cgx_web[1], cgy_web[1], cgz_web[1], color=colors[1], alpha=1, marker='.', linewidths=0, s=1);
ax.plot3D(cgx_web[1], cgy_web[1], cgz_web[1], color = colors[1], linestyle = ":", linewidth=0.9, label = "{}".format(webs[1]))
#sc = ax.scatter(cgx_sat_all[0], cgy_sat_all[0], cgz_sat_all[0], color=colors[2], alpha=1, marker='+', linewidths=0, s=1);
ax.plot3D(cgx_sat_all[0], cgy_sat_all[0], cgz_sat_all[0], color = colors[2], linestyle = ":", linewidth=0.9, label = "{}".format(sats[0]))
#sc = ax.scatter(cgx_sat_all[1], cgy_sat_all[1], cgz_sat_all[1], color=colors[3], alpha=1, marker='+', linewidths=0, s=1);
ax.plot3D(cgx_sat_all[1], cgy_sat_all[1], cgz_sat_all[1], color = colors[3], linestyle = ":", linewidth=0.9, label = "{}".format(sats[1]))
for i in np.linspace(0, len(redshift_array)-1, 6).astype(int):
    ax.text(cgx_sat_all[1][i]-0.0075,cgy_sat_all[1][i]-0.0075,cgz_sat_all[1][i]+0.015,'z={}'.format(round(redshift_array[i],2)), size=4, zorder=1,color='k')
    ax.text(cgx_sat_all[0][i]+0.003,cgy_sat_all[0][i]+0.003,cgz_sat_all[0][i]-0.015,'z={}'.format(round(redshift_array[i],2)), size=4, zorder=1,color='k')         
plt.xlim(0.225,0.38)
plt.ylim(0.425,0.58)
ax.set_zlim3d(0.175,0.33)
ax.set_ylabel(r'$y$', fontsize=10)
ax.set_zlabel(r'$z$', fontsize=10)
ax.tick_params(axis="x", labelsize=8)
ax.tick_params(axis="y", labelsize=8)
ax.tick_params(axis="z", labelsize=8)
plt.legend(loc = "upper left",fontsize=4)
plt.title('Simulación para centros de masa')
ax.view_init(0, 0)
plt.savefig(images_path + 'Simulación para centros de masa, 0, 0.png', dpi=500)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#sc = ax.scatter(cgx_web[0], cgy_web[0], cgz_web[0], color=colors[0], alpha=1, marker='.', linewidths=0, s=1);
ax.plot3D(cgx_web[0], cgy_web[0], cgz_web[0], color = colors[0], linestyle = ":", linewidth=0.9, label = "{}".format(webs[0]))
#sc = ax.scatter(cgx_web[1], cgy_web[1], cgz_web[1], color=colors[1], alpha=1, marker='.', linewidths=0, s=1);
ax.plot3D(cgx_web[1], cgy_web[1], cgz_web[1], color = colors[1], linestyle = ":", linewidth=0.9, label = "{}".format(webs[1]))
#sc = ax.scatter(cgx_sat_all[0], cgy_sat_all[0], cgz_sat_all[0], color=colors[2], alpha=1, marker='+', linewidths=0, s=1);
ax.plot3D(cgx_sat_all[0], cgy_sat_all[0], cgz_sat_all[0], color = colors[2], linestyle = ":", linewidth=0.9, label = "{}".format(sats[0]))
#sc = ax.scatter(cgx_sat_all[1], cgy_sat_all[1], cgz_sat_all[1], color=colors[3], alpha=1, marker='+', linewidths=0, s=1);
ax.plot3D(cgx_sat_all[1], cgy_sat_all[1], cgz_sat_all[1], color = colors[3], linestyle = ":", linewidth=0.9, label = "{}".format(sats[1]))
for i in np.linspace(0, len(redshift_array)-1, 6).astype(int):
    ax.text(cgx_sat_all[1][i]-0.015,cgy_sat_all[1][i]-0.015,cgz_sat_all[1][i],'z={}'.format(round(redshift_array[i],2)), size=4, zorder=1,color='k')
    ax.text(cgx_sat_all[0][i]+0.015,cgy_sat_all[0][i]+0.015,cgz_sat_all[0][i],'z={}'.format(round(redshift_array[i],2)), size=4, zorder=1,color='k')         
plt.xlim(0.225,0.38)
plt.ylim(0.425,0.58)
ax.set_zlim3d(0.175,0.33)
ax.set_xlabel(r'$x$', fontsize=10)
ax.set_ylabel(r'$y$', fontsize=10)
ax.tick_params(axis="x", labelsize=8)
ax.tick_params(axis="y", labelsize=8)
ax.tick_params(axis="z", labelsize=8)
plt.legend(loc = "upper left",fontsize=4)
plt.title('Simulación para centros de masa')
ax.view_init(90, 90)
plt.savefig(images_path + 'Simulación para centros de masa, 90, 90.png', dpi=500)
plt.show()

nweb2_array = []

for n in np.linspace(0,60,61).astype(int):
    nweb2 = n_web_array[int(n)]
    nweb2_array.append(nweb2)
    
nweb3_array = []

for n in np.linspace(60,121,61).astype(int):
    nweb3 = n_web_array[int(n)]
    nweb3_array.append(nweb3)

fig = plt.figure()
plt.title(r'Número de partículas corregido')
plt.xlabel('$t$')
plt.ylabel('$n_T^*(t)$')
plt.axis([min(atime_array), max(atime_array), 23000, 42000])
plt.axhline(y=n_tot_web_array[0], xmin=0, xmax=max(atime_array), color='b', linewidth = 1)
plt.axhline(y=n_tot_web_array[1], xmin=0, xmax=max(atime_array), color='r', linewidth = 1)
plt.plot(atime_array, nweb2_array, color = "b", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, nweb2_array, color='b', marker='+', s=5e1, linewidths=1, label = r'$n_{T,web2}^*$')
plt.plot(atime_array, nweb3_array, color = "r", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, nweb3_array, color='r', marker='+', s=5e1, linewidths=1, label = r'$n_{T,web3}^*$')
plt.legend(loc='lower center')
plt.grid(True, which="both", ls="-")
plt.savefig(images_path + 'Evolución del número corregido.png', dpi=700)
plt.show()

end = time.time()
print('___________________________________________________________________________________________________________')
print("Tiempo de ejecución (s):", end-start)
print("En minutos (min):", (end-start)/60)
print('___________________________________________________________________________________________________________')

