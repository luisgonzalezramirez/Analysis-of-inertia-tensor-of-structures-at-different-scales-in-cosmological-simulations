
#Código de análisis de la convergencia de centros de masa:
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

frames = listdir(directory)

webs = listdir(structure)

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
dist_meanarray = [] #Array para distancias medias entre centro de masas y partículas.

cgx_array = []
cgy_array = []
cgz_array = []

cgx_array_con = []
cgy_array_con = []
cgz_array_con = []
cg_array_con = []

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
    
    print('Redshift, z =',z)
    
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
    
    #Se definen arrays para las estructuras web.
    
    rmweb_array = [] #Masas de estructuras web.
    xweb_array = [] #Posición x estructura web.
    yweb_array = [] #Posición y estructura web.
    zweb_array = [] #Posición z estructura web.
    
    n = 0
    
    end = time.time()
    print('___________________________________________________________________________________________________________')
    print("Tiempo de ejecución en punto 2 (s):", end-start)
    print("En minutos en punto 2 (min):", (end-start)/60)
    print('___________________________________________________________________________________________________________')

    dist_mean_array = []
    dist_meanarray = []
    
    #Se obtienen los arrays de las estructuras web a partir de array de dos arrays con los números de partículas:
    for web in web_arrayconverted:
        for j in web:
            mweb = m_array[int(j)-1]
            rmweb_array.append(mweb)
            xweb = x_array[int(j)-1]
            xweb_array.append(xweb)
            yweb = y_array[int(j)-1]
            yweb_array.append(yweb)
            zweb = z_array[int(j)-1]
            zweb_array.append(zweb)
        
        #Se calculan los centros de masa iniciales:
        cgx_i = np.sum(np.multiply(xweb_array,rmweb_array))/np.sum(rmweb_array)
        cgy_i = np.sum(np.multiply(yweb_array,rmweb_array))/np.sum(rmweb_array)
        cgz_i = np.sum(np.multiply(zweb_array,rmweb_array))/np.sum(rmweb_array)
        cg_vector = [webs[n],z,cgx_i,cgy_i,cgz_i]
        cg_array.append(cg_vector)
        n = n + 1
        
    dist_array = []
    dist_generalarray = []
    dist_max_array = []
    dist_generalmaxarray = []

    n = 0
    
    #Se calcula las distancias de las partículas de cada estructura al centro de masas inicial: 
    for web in web_arrayconverted:
        dist_mean_array = [webs[n],z]
        dist_max_array = [webs[n],z]
        for j in web:
            xweb = x_array[int(j)-1]
            yweb = y_array[int(j)-1]
            zweb = z_array[int(j)-1]
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

#Se grafica en una misma imagen la trayectoria de los centros de masa en convergencia:
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for m in np.linspace(0, len(cg_array)-1, len(cg_array)).astype(int):
    if cg_array[m][0] == webs[0]:
        ax.scatter(cg_array[m][2], cg_array[m][3], cg_array[m][4], color='r', marker='+', s=1e2,linewidths=0.25)
    elif cg_array[m][0] == webs[1]:
        ax.scatter(cg_array[m][2], cg_array[m][3], cg_array[m][4], color='b', marker='+', s=1e2,linewidths=0.25)
plt.xlim(0.35, 0.45)
plt.ylim(0.25, 0.35)
ax.set_zlim3d(0.35, 0.45)
plt.title('Zoom Simulación convergencia')
plt.savefig(images_path + 'Zoom Simulación.png', dpi=500)
plt.close() 
    
f.close()
    
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
    redshift_z = (1-atime)/atime
    
    print('Redshift, z =',redshift_z)
    
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
    
    n = 0
    
    end = time.time()
    print('___________________________________________________________________________________________________________')
    print("Tiempo de ejecución en punto 2 (s):", end-start)
    print("En minutos en punto 2 (min):", (end-start)/60)
    print('___________________________________________________________________________________________________________')

    i = 0
    k = 0
    n = 0

    #Se obtienen los arrays de las estructuras web a partir de array de dos arrays con los números de partículas:
    for web in web_arrayconverted:
        #Se definen arrays para las estructuras web.
        rmweb_array_con = [] #Masas de estructuras web.
        xweb_array_con = [] #Posición x estructura web.
        yweb_array_con = [] #Posición y estructura web.
        zweb_array_con = [] #Posición z estructura web.
        for j in web:
            mweb_con = m_array[int(j)-1]
            rmweb_array_con.append(mweb_con)
            xweb_con = x_array[int(j)-1]
            xweb_array_con.append(xweb_con)
            yweb_con = y_array[int(j)-1]
            yweb_array_con.append(yweb_con)
            zweb_con = z_array[int(j)-1]
            zweb_array_con.append(zweb_con)
        cgx = np.sum(np.multiply(xweb_array_con,rmweb_array_con))/np.sum(rmweb_array_con)
        cgy = np.sum(np.multiply(yweb_array_con,rmweb_array_con))/np.sum(rmweb_array_con)
        cgz = np.sum(np.multiply(zweb_array_con,rmweb_array_con))/np.sum(rmweb_array_con)
        xweb_array_new = xweb_array_con
        yweb_array_new = yweb_array_con
        zweb_array_new = zweb_array_con
        k = 0
        for (x,y,z) in zip(xweb_array_con,yweb_array_con,zweb_array_con):
            p1 = np.array([x, y, z])
            p2 = np.array([cgx, cgy, cgz])
            squared_dist = np.sum((p1-p2)**2, axis=0)
            dist = np.sqrt(squared_dist)
            #Se determina el límite para que se sitúen las partículas en el 95% del radio máximo:
            if dist >= dist_generalarray[i][2]/(dist_generalmaxarray[i][2]*0.95):
                del xweb_array_new[k]
                del yweb_array_new[k]
                del zweb_array_new[k]
                del rmweb_array_con[k]
            k = k + 1
        cgx_con = np.sum(np.multiply(xweb_array_new,rmweb_array_con))/np.sum(rmweb_array_con)
        cgy_con = np.sum(np.multiply(yweb_array_new,rmweb_array_con))/np.sum(rmweb_array_con)
        cgz_con = np.sum(np.multiply(zweb_array_new,rmweb_array_con))/np.sum(rmweb_array_con)
        
        i = i + 1
        
        cg_vector_con = [webs[n],redshift_z,cgx_con,cgy_con,cgz_con]
        cg_array_con.append(cg_vector_con)
        
        cgx_array_con.append(cgx_con)
        cgy_array_con.append(cgy_con)
        cgz_array_con.append(cgz_con)
        
        n = n + 1

#Se grafica en una misma imagen la trayectoria de los centros de masa en convergencia:
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for m in np.linspace(0, len(cg_array)-1, len(cg_array)).astype(int):
    if cg_array_con[m][0] == webs[0]:
        ax.scatter(cg_array_con[m][2], cg_array_con[m][3], cg_array_con[m][4], color='r', marker='+', s=1e2,linewidths=0.25)
    elif cg_array[m][0] == webs[1]:
        ax.scatter(cg_array_con[m][2], cg_array_con[m][3], cg_array_con[m][4], color='b', marker='+', s=1e2,linewidths=0.25)
plt.xlim(0.35, 0.45)
plt.ylim(0.25, 0.35)
ax.set_zlim3d(0.35, 0.45)
plt.title('Zoom Simulación convergencia z = {:.2f}'.format(redshift_z))
plt.savefig(images_path + 'Zoom Simulación, con convergencia.png', dpi=500)
plt.close()        
        
end = time.time()
print('___________________________________________________________________________________________________________')
print("Tiempo de ejecución final (s):", end-start)
print("En minutos (min):", (end-start)/60)
print('___________________________________________________________________________________________________________')