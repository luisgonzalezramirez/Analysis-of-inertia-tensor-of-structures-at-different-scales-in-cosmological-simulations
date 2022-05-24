
#Función para el cálculo de las partículas de los satélites:
#___________________________________________________________#

import os
import shutil
import numpy as np
import time
from os import listdir

start = time.time()

#Input de los directorios para separar los archivos:
print('___________________________________________________________________________________________________________')
target = input('Introduce dirección local del código:\n')

print('___________________________________________________________________________________________________________')
satellites_path = input('Introduce dirección local de entrada de los archivos de satélites:\n')

print('___________________________________________________________________________________________________________')
satellites_simulation_path = input('Introduce dirección local del archivo d5004.41132:\n')

for root, directories, sats in os.walk(satellites_path, topdown=False):
    for name in sats:
        f = os.path.join(root, name)
        shutil.copy2(f, target)

for root, directories, sats_sim in os.walk(satellites_simulation_path, topdown=False):
    for name in sats_sim:
        f = os.path.join(root, name)
        shutil.copy2(f, target)

frame_sat = f

n = 0

sat_data = []

#Obtención de los centros y radios de los satélites:
for sat in sats:
    s = open(sat,'r')

    line_count = 0
    for l in s:
        if l != "\n":
            line_count += 1
            
    length = line_count
    
    s.close()
    
    s = open(sat,'r')
    
    r = range(0,line_count,1)
    
    sat_data.append(sats[n])
    
    #Se definen los arrays de los centros de masa de los satélites en coordenadas x, y, z, y el radio r:
    cx_array = []
    cy_array = []
    cz_array = []
    r_array = []
    sat_array = []
    
    #Se lee cada línea (líneas de 4 caracteres correspondientes al centro, 3 primeras, y radio) y se obtiene un array de
    #dos arrays con las satélites persistentes y no persistentes:
    for line in r:
        data_line = s.readline()
        cx,cy,cz,r = data_line.split()
        cx_array.append(cx)
        cy_array.append(cy)
        cz_array.append(cz)
        r_array.append(r)
        
	#Tarda más ponerlo en una matriz, pero visualmente es útil o si se quiere exportar un Excel:
    sat_array = np.array([cx_array,cy_array,cz_array,r_array])
    sat_data.append(sat_array)
    n = n + 1

#Aquí se pueden leer y se calcular:
#Satelites(v_1x,v_1y,v_1z,v_2x,v_2y,v_2z,v_data,sat_directory,target_directory)

sats_directory = satellites_path
target_directory = target

import numpy as np
import struct

sats = listdir(sats_directory)

f = open(frame_sat, "rb")



#for root, directories, sats in os.walk(target_directory, topdown=False):
#    for name in sats:
#        f = os.path.join(root, name)
#        shutil.copy2(f, target_directory)

#Lectura de parámetros:
parameters1 = f.read(20)
l1, out, itime, irun, l2 = struct.unpack(">iiiii", parameters1)

parameters2 = f.read(20)
l1, out, itstop, nobj, l2 = struct.unpack(">iiiii", parameters2)


parameters3 = f.read(20)
l1, out, itdump, ngas, l2 = struct.unpack(">iiiii", parameters3)

parameters4 = f.read(20)
l1, out, iout, ndark, l2 = struct.unpack(">iiiii", parameters4)

parameters5 = f.read(20)
l1, out, rtime, L, l2 = struct.unpack(">iiffi", parameters5)

parameters6 = f.read(20)
l1, out, atime, intl, l2 = struct.unpack(">iifii", parameters6)

#Redshift:
z = (1-atime)/atime

print('Redshift, z =',z)

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


v_1x = x_array
v_1y = y_array
v_1z = z_array
v_type = type_array

v_data = sat_data

def Satélites(v_1x,v_1y,v_1z,v_type,v_2x,v_2y,v_2z,v_data,sats_directory,target_directory,frame_sat):
    
    '''
    Función:
    ________________________________________________________________________________________
    Realiza la lectura de los datos de la simulación d5004.41132 para obtener las partículas de los satélites.
    
    Variables:
    ________________________________________________________________________________________
    * v_1x, variable correspondiente al array de las posiciones x de las partículas (x_array).
    * v_1y, variable correspondiente al array de las posiciones y de las partículas (y_array).
    * v_1z, variable correspondiente al array de las posiciones z de las partículas (z_array).
    * v_type, variable correspondiente al array de los tipos de las partículas.
    * v_2x, variable del array de las coordenadas x de los centros de satélites (cx_array).
    * v_2y, variable del array de las coordenadas y de los centros de satélites (cy_array).
    * v_2z, variable del array de las coordenadas z de los centros de satélites (cz_array).
    * v_data, array de los datos de los satélites (sat_data).
    * sats_directory, directorio de los satélites (satellites_path).
    * target_directory, directorio de del código (target).
    * frame_sat, 
    '''
    
    rbis = range(0,len(v_1x),1)

    i = 0
    m = 0

    for n in [1, 3]:
        cx_array = v_data[n][0,:]
        cy_array = v_data[n][1,:]
        cz_array = v_data[n][2,:]
        r_array = v_data[n][3,:]
        for i in np.linspace(0, len(r_array)-1, len(r_array)).astype(int):
            center = np.array([cx_array[i], cy_array[i], cz_array[i]])
            center = center.astype(np.float)
            par_sat_array_info = [v_data[n-1],i]
            par_sat_array = []
            for j in rbis:
                x = v_1x[int(j)]
                y = v_1y[int(j)]
                z = v_1z[int(j)]
                typ = v_type[int(j)]
                p = np.array([x, y, z])
                squared_dist_sat = np.sum((center-p)**2, axis=0)
                dist_sat = np.sqrt(squared_dist_sat)
                if dist_sat <= float(r_array[i]):
                    if typ != 0:
                        par_sat_array_info.append(j)
                        par_sat_array.append(j)
            f_info = open('{}{}.txt'.format(v_data[n-1],i),'w+')
            f = open('Sat, {}{}.txt'.format(v_data[n-1],i),'w+')
            for m in par_sat_array:
                f_info.write('{}\n'.format(m))
                f.write('{}\n'.format(m))
            f.close()

Satélites(v_1x,v_1y,v_1z,v_type,cx_array,cy_array,cz_array,sat_data,sats_directory,target_directory,frame_sat)


end = time.time()
print('___________________________________________________________________________________________________________')
print("Tiempo de ejecución (s):", end-start)
print("En minutos (min):", (end-start)/60)
print('___________________________________________________________________________________________________________')