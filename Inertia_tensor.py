
#Código de lectura de datos y cálculo del tensor de inercia:
#_______________________________________________________________#

#Importar módulos de libreria de python:

import time
import os
from os import listdir
import shutil
import struct
import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
import statistics as st

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
#Se define una lista para los eigenvectores:
evecs_array_array = []

m = 0

num_i = 0 

T_web_array_array = []
e_web_array_array = []
p_web_array_array = []

c_b_array_array = []
c_c_array_array = []

web_length_array = []

#Obtención de partículas de las estructuras web2 y web3, correspondientes a un filamento y un plano:
for web in webs:
    
    T_web_array_array.append(webs[m])
    e_web_array_array.append(webs[m])
    p_web_array_array.append(webs[m])
    
    c_b_array_array.append(webs[m])
    c_c_array_array.append(webs[m])
    
    #Se definen arrays para los semiejes:
    a_array = []
    b_array = []
    c_array = []
    
    #Se definen arrays para los parámetros:
    T_web_array = []
    e_web_array = []
    p_web_array = []
    
    evecs_array = [webs[m]]
    
    c_b_array = []
    c_c_array = []
    
    #a_array.append(webs[m])
    #b_array.append(webs[m])
    #c_array.append(webs[m])
    
    fweb = open(web, "r")
    #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
    num_par = fweb.readlines()
    web_length = num_par[0]
    web_length_array.append(web_length)
    #Se ha definido una lista nueva para el número de partículas eliminando la posición 0 que es el valor de la longitud:
    del num_par[0]
    length = len(num_par)
    rweb = range(len(num_par))

    #Se definen arrays:
    array_data = [] #Array para la posibilidad de acumular todo.
    cg_array = [] #Centro de masa.
    dist_generalarray = [] #Array para distancias entre centro de masas y partículas.

    cgx_array = []
    cgy_array = []
    cgz_array = []
    
    redshift_array = []
    atime_array = []
    
    
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
        rmweb_array = [] #Masas de estructuras web.
        xweb_array = [] #Posición x estructura web.
        yweb_array = [] #Posición y estructura web.
        zweb_array = [] #Posición z estructura web.

        k = 0

        #Se obtienen los arrays de las estructuras web a partir de array de dos arrays con los números de partículas:
        for j in num_par:
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
        cg_vector = [webs[k],z,cgx_i,cgy_i,cgz_i]
        cg_array.append(cg_vector)
        k = k + 1
            
        dist_array = []
        dist_generalarray = []
        dist_max_array = []
        dist_generalmaxarray = []

        #Se calcula las distancias de las partículas de cada estructura al centro de masas inicial: 
        dist_mean_array = [webs[m],z]
        dist_max_array = [webs[m],z]
        for j in num_par:
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

        rmweb_array_con = [] #Masas de estructuras web.
        xweb_array_con = [] #Posición x estructura web.
        yweb_array_con = [] #Posición y estructura web.
        zweb_array_con = [] #Posición z estructura web.

        #Se obtienen los arrays de las estructuras web a partir de array de dos arrays con los números de partículas:
        for j in num_par:
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
            if dist >= dist_mean_array[2]/(dist_max_array[2]*0.95):
                del xweb_array_new[k]
                del yweb_array_new[k]
                del zweb_array_new[k]
                del rmweb_array_con[k]
            k = k + 1
        cgx = np.sum(np.multiply(xweb_array_new,rmweb_array_con))/np.sum(rmweb_array_con)
        cgy = np.sum(np.multiply(yweb_array_new,rmweb_array_con))/np.sum(rmweb_array_con)
        cgz = np.sum(np.multiply(zweb_array_new,rmweb_array_con))/np.sum(rmweb_array_con)

        cg_vector_con = [webs[m],z,cgx,cgy,cgz]
        
        cgx_array.append(cgx)
        cgy_array.append(cgy)
        cgz_array.append(cgz)
        
        evals,evecs,a,b,c = Componentes_Web_no_print(cgx,cgy,cgz,xweb_array,yweb_array,zweb_array,rmweb_array)
        
        M_total_web = np.sum(evals)/2
        
        M_total_web_strc = np.sum(rmweb_array)
        
        evecs_array.append(evecs)
        
        c_b = b/a
        c_c = c/a
        
        c_b_array.append(c_b)
        c_c_array.append(c_c)
        
        a_array.append(a)
        b_array.append(b)
        c_array.append(c)
        
        #Parámetro de triaxialidad:
        T_web = (1 - (b**2)/(a**2))/(1 - (c**2)/(a**2))
        T_web_array.append(T_web)
        
        #Elipticidad:
        e_web = (a**2 - c**2)/(a**2 + b**2 + c**2)
        e_web_array.append(e_web)
        
        #Prolaticidad:
        p_web = (a**2 + c**2 - 2*b**2)/(a**2 + b**2 + c**2)
        p_web_array.append(p_web)
    
    num_i = num_i + 1 
    
    c_b_array_array.append(c_b_array)
    c_c_array_array.append(c_c_array)
    
    T_web_array_array.append(T_web_array)
    e_web_array_array.append(e_web_array)
    p_web_array_array.append(p_web_array)
    
    evecs_array_array.append(evecs_array)
    
    figp1 = plt.figure()
    plt.title('Triaxialidad: {}'.format(webs[m]))
    plt.ylabel('T(t)')
    plt.axis([min(atime_array), max(atime_array), min(T_web_array)-min(T_web_array)*0.1, 1])
    plt.plot(atime_array, T_web_array, color = "red", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, T_web_array, color='r', marker='+', s=1e2, linewidths=1, label = "T(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Triaxialidad {}.png'.format(webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    figp2 = plt.figure()
    plt.title('Elipticidad: {}'.format(webs[m]))
    plt.ylabel('e(t)')
    plt.axis([min(atime_array), max(atime_array), min(e_web_array)-min(e_web_array)*0.1, 1])
    plt.plot(atime_array, e_web_array, color = "blue", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, e_web_array, color='b', marker='+', s=1e2, linewidths=1, label = "e(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Elipticidad {}.png'.format(webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    figp3 = plt.figure()
    plt.title('Prolaticidad: {}'.format(webs[m]))
    plt.ylabel('p(t)')
    plt.axis([min(atime_array), max(atime_array), min(p_web_array)-min(p_web_array)*0.1, 1])
    plt.plot(atime_array, p_web_array, color = "green", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, p_web_array, color='g', marker='+', s=1e2, linewidths=1, label = "p(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Prolaticidad {}.png'.format(webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig = plt.figure()
    plt.title('Estructura: {}'.format(webs[m]))
    plt.ylabel('Semiejes')
    plt.axis([max(redshift_array), min(redshift_array), 0, max(a_array)+max(a_array)*0.1])
    plt.plot(redshift_array, a_array, color = "red", linewidth = 1, linestyle = "-", label = "a(z)")
    plt.plot(redshift_array, b_array, color = "blue", linewidth = 1, linestyle = "-", label = "b(z)")
    plt.plot(redshift_array, c_array, color = "green", linewidth = 1, linestyle = "-", label = "c(z)")
    plt.xlabel('z')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Semiejes {}.png'.format(webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig2 = plt.figure()
    plt.title('Estructura: {}'.format(webs[m]))
    plt.ylabel('Semiejes')
    plt.axis([min(atime_array), max(atime_array), 0, max(a_array)+max(a_array)*0.1])
    plt.plot(atime_array, a_array, color = "red", linewidth = 1, linestyle = "-", label = "a(t)")
    plt.plot(atime_array, b_array, color = "blue", linewidth = 1, linestyle = "-", label = "b(t)")
    plt.plot(atime_array, c_array, color = "green", linewidth = 1, linestyle = "-", label = "c(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Semiejes, atime, {}.png'.format(webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig3 = plt.figure()
    plt.title('Estructura: {}'.format(webs[m]))
    plt.ylabel('Cocientes de semiejes')
    plt.axis([min(atime_array), max(atime_array), 0, max(c_b_array)+max(c_c_array)*0.1])
    plt.plot(atime_array, c_b_array, color = "red", linewidth = 1, linestyle = ":")
    plt.plot(atime_array, c_c_array, color = "blue", linewidth = 1, linestyle = ":")
    
    line_t = np.linspace(min(atime_array), max(atime_array), 100)
    
    cb_fit = np.poly1d(np.polyfit(atime_array, c_b_array, 5))
    cc_fit = np.poly1d(np.polyfit(atime_array, c_c_array, 5))

    plt.scatter(atime_array, c_b_array, color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t)")
    plt.scatter(atime_array, c_c_array, color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t)")
    plt.plot(line_t, cb_fit(line_t),color = 'r')
    plt.plot(line_t, cc_fit(line_t),color = 'b')
    
    plt.grid()
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.savefig(images_path + 'Semiejes, atime, cocientes, {}.png'.format(webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig4 = plt.figure()
    plt.title('Estructura: {}'.format(webs[m]))
    plt.ylabel('Cocientes de semiejes')
    plt.axis([max(redshift_array), min(redshift_array), 0, max(c_b_array)+max(c_c_array)*0.1])
    plt.plot(redshift_array, c_b_array, color = "red", linewidth = 1, linestyle = ":")
    plt.plot(redshift_array, c_c_array, color = "blue", linewidth = 1, linestyle = ":")
    print(webs[m],'A z = 9, b/a =',c_b_array[0])
    print(webs[m],'A z = 9, c/a =',c_c_array[0])
    print(webs[m],'A z = 4, b/a =',c_b_array[3])
    print(webs[m],'A z = 4, c/a =',c_c_array[3])
    print(webs[m],'A z = 2, b/a =',c_b_array[14])
    print(webs[m],'A z = 2, c/a =',c_c_array[14])
    print(webs[m],'A z = 0, b/a =',c_b_array[60])
    print(webs[m],'A z = 0, c/a =',c_c_array[60])
    plt.scatter(redshift_array, c_b_array, color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t)")
    plt.scatter(redshift_array, c_c_array, color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t)")
    
    plt.grid()
    plt.xlabel('z')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.savefig(images_path + 'Semiejes, atime, cocientes, {}.png'.format(webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    m = m + 1

pickle.dump(evecs_array_array, open("evecs_array_array.dat", "wb"))
    
#Se define una lista para los redshift:
redshift_array = []
#Se define una lista para los eigenvectores:
evecs_array_array_sat = []

m = 0

T_sat_array_array = []
e_sat_array_array = []
p_sat_array_array = []

c_b_array_array_sat = []
c_c_array_array_sat = []

sat_length_array = []

k = 0

#Obtención de partículas de los satélites persistentes y no persistentes:
for sat in sats:
    
    T_sat_array_array.append(sats[m])
    e_sat_array_array.append(sats[m])
    p_sat_array_array.append(sats[m])
    
    c_b_array_array_sat.append(sats[m])
    c_c_array_array_sat.append(sats[m])
    
    #Se definen arrays para los semiejes:
    a_array_sat = []
    b_array_sat = []
    c_array_sat = []
    
    #Se definen arrays para los parámetros:
    T_sat_array = []
    e_sat_array = []
    p_sat_array = []
    
    evecs_array_sat = [sats[m]]
    
    c_b_array_sat = []
    c_c_array_sat = []
    
    #Se definen arrays:
    array_data_sat = [] #Array para la posibilidad de acumular todo.
    cg_array_sat = [] #Centro de masa.
    
    cgx_array_sat = []
    cgy_array_sat = []
    cgz_array_sat = []
    
    redshift_array = []
    
    sat_arrayconverted = []
    sat_array = []
    
    #Obtención de partículas de los satélites, persistentes y no persistentes:
    fsat = open(sat, "r")
    #Se lee cada línea (líneas de un caracter correspondiente al número de partícula):
    num_par_sat = fsat.readlines()
    sat_length = len(num_par_sat)
    sat_length_array.append(sat_length)
    #Se crean un array para los números de las partículas de cada estructura:
    converted_num_par_sat_array = np.array(num_par_sat)
    converted_num_par_sat_array = converted_num_par_sat_array.astype(int)
    #Se obtiene un array de dos arrays con las partículas de las dos grupos de satélites:
    sat_arrayconverted.append(converted_num_par_sat_array)
    sat_array.append(num_par)
    
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
        cg_vector_sat = [sat[k],z,cgx_sat,cgy_sat,cgz_sat]
        cg_array_sat.append(cg_vector_sat)
        
        evals_sat,evecs_sat,a_sat,b_sat,c_sat = Componentes_Web_no_print(cgx_sat,cgy_sat,cgz_sat,xsat_array,ysat_array,zsat_array,rmsat_array)
        
        M_total = np.sum(evals_sat)/2
        
        M_total_strc = np.sum(rmsat_array)
        
        evecs_array_sat.append(evecs_sat)
        
        c_b_sat = b_sat/a_sat
        c_c_sat = c_sat/a_sat
        
        c_b_array_sat.append(c_b_sat)
        c_c_array_sat.append(c_c_sat)
        
        a_array_sat.append(a_sat)
        b_array_sat.append(b_sat)
        c_array_sat.append(c_sat)
        
        #Parámetro de triaxialidad:
        T_sat = (1 - (b_sat**2)/(a_sat**2))/(1 - (c_sat**2)/(a**2))
        T_sat_array.append(T_sat)
        
        #Elipticidad:
        e_sat = (a_sat**2 - c_sat**2)/(a_sat**2 + b_sat**2 + c_sat**2)
        e_sat_array.append(e_sat)
        
        #Prolaticidad:
        p_sat = (a_sat**2 + c_sat**2 - 2*b_sat**2)/(a_sat**2 + b_sat**2 + c_sat**2)
        p_sat_array.append(p_sat)
    
    k = k + 1
    
    c_b_array_array_sat.append(c_b_array_sat)
    c_c_array_array_sat.append(c_c_array_sat)
    
    T_sat_array_array.append(T_sat_array)
    e_sat_array_array.append(e_sat_array)
    p_sat_array_array.append(p_sat_array)
    
    evecs_array_array_sat.append(evecs_array_sat)
    
    figp1 = plt.figure()
    plt.title('Triaxialidad: {}'.format(sats[m]))
    plt.ylabel('T(t)')
    plt.axis([min(atime_array), max(atime_array), min(T_sat_array)-min(T_sat_array)*0.1, 1])
    plt.plot(atime_array, T_sat_array, color = "red", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, T_sat_array, color='r', marker='+', s=1e2, linewidths=1, label = "T(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Triaxialidad {}.png'.format(sats[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    figp2 = plt.figure()
    plt.title('Elipticidad: {}'.format(sats[m]))
    plt.ylabel('e(t)')
    plt.axis([min(atime_array), max(atime_array), min(e_sat_array)-min(e_sat_array)*0.1, 1])
    plt.plot(atime_array, e_sat_array, color = "blue", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, e_sat_array, color='b', marker='+', s=1e2, linewidths=1, label = "e(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Elipticidad {}.png'.format(sats[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    figp3 = plt.figure()
    plt.title('Prolaticidad: {}'.format(sats[m]))
    plt.ylabel('p(t)')
    plt.axis([min(atime_array), max(atime_array), min(p_sat_array)-min(p_sat_array)*0.1, 1])
    plt.plot(atime_array, p_sat_array, color = "green", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, p_sat_array, color='g', marker='+', s=1e2, linewidths=1, label = "p(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Prolaticidad {}.png'.format(sats[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig = plt.figure()
    plt.title('Estructura: {}'.format(sats[m]))
    plt.ylabel('Semiejes')
    plt.axis([max(redshift_array), min(redshift_array), 0, max(a_array_sat)+max(a_array_sat)*0.1])
    plt.plot(redshift_array, a_array_sat, color = "red", linewidth = 1, linestyle = "-", label = "a(z)")
    plt.plot(redshift_array, b_array_sat, color = "blue", linewidth = 1, linestyle = "-", label = "b(z)")
    plt.plot(redshift_array, c_array_sat, color = "green", linewidth = 1, linestyle = "-", label = "c(z)")
    plt.xlabel('z')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Semiejes {}.png'.format(sats[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig2 = plt.figure()
    plt.title('Estructura: {}'.format(sats[m]))
    plt.ylabel('Semiejes')
    plt.axis([min(atime_array), max(atime_array), 0, max(a_array_sat)+max(a_array_sat)*0.1])
    plt.plot(atime_array, a_array_sat, color = "red", linewidth = 1, linestyle = "-", label = "a(t)")
    plt.plot(atime_array, b_array_sat, color = "blue", linewidth = 1, linestyle = "-", label = "b(t)")
    plt.plot(atime_array, c_array_sat, color = "green", linewidth = 1, linestyle = "-", label = "c(t)")
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Semiejes, atime, {}.png'.format(sats[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig3 = plt.figure()
    plt.title('Estructura: {}'.format(sats[m]))
    plt.ylabel('Cocientes de semiejes')
    plt.axis([min(atime_array), max(atime_array), 0, max(c_b_array_sat)+max(c_c_array_sat)*0.1])
    plt.plot(atime_array, c_b_array_sat, color = "red", linewidth = 1, linestyle = ":")
    plt.plot(atime_array, c_c_array_sat, color = "blue", linewidth = 1, linestyle = ":")
    
    
    line_t = np.linspace(min(atime_array), max(atime_array), 100)
    
    cb_fit_sat = np.poly1d(np.polyfit(atime_array, c_b_array_sat, 5))
    cc_fit_sat = np.poly1d(np.polyfit(atime_array, c_c_array_sat, 5))

    plt.scatter(atime_array, c_b_array_sat, color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t)")
    plt.scatter(atime_array, c_c_array_sat, color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t)")
    plt.plot(line_t, cb_fit_sat(line_t),color = 'r')
    plt.plot(line_t, cc_fit_sat(line_t),color = 'b')
    
    plt.grid()
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.savefig(images_path + 'Semiejes, atime, cocientes, {}.png'.format(sats[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    fig4 = plt.figure()
    plt.title('Estructura: {}'.format(sats[m]))
    plt.ylabel('Cocientes de semiejes')
    plt.axis([max(redshift_array), min(redshift_array), 0, max(c_b_array_sat)+max(c_c_array_sat)*0.1])
    plt.plot(redshift_array, c_b_array_sat, color = "red", linewidth = 1, linestyle = ":")
    plt.plot(redshift_array, c_c_array_sat, color = "blue", linewidth = 1, linestyle = ":")

    plt.scatter(redshift_array, c_b_array_sat, color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t)")
    plt.scatter(redshift_array, c_c_array_sat, color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t)")
    print(sats[m],'A z = 9, b/a =',c_b_array_sat[0])
    print(sats[m],'A z = 9, b/a =',c_c_array_sat[0])
    print(sats[m],'A z = 4, b/a =',c_b_array_sat[3])
    print(sats[m],'A z = 4, c/a =',c_c_array_sat[3])
    print(sats[m],'A z = 2, b/a =',c_b_array_sat[14])
    print(sats[m],'A z = 2, b/a =',c_c_array_sat[14])
    print(sats[m],'A z = 0, b/a =',c_b_array_sat[60])
    print(sats[m],'A z = 0, b/a =',c_c_array_sat[60])
    plt.grid()
    plt.xlabel('z')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.savefig(images_path + 'Semiejes, z, cocientes, {}.png'.format(sats[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    m = m + 1


pickle.dump(evecs_array_array_sat, open("evecs_array_array_sat.dat", "wb"))

#Aquí se pueden guardar previamente los arrays de los autovectores:

m = 0

pickle.dump(atime_array, open("atime_array.dat", "wb"))
pickle.dump(atime_array, open("redshift_array.dat", "wb"))

evecs_array_array = pickle.load(open("evecs_array_array.dat", "rb"))
evecs_array_array_sat = pickle.load(open("evecs_array_array_sat.dat", "rb"))

evecs_array_array_sat_element1 = evecs_array_array_sat[1]
evecs_array_array_sat_element2 = evecs_array_array_sat[0]

evecs_array_array_sat_order = [evecs_array_array_sat_element1, evecs_array_array_sat_element2]
evecs_array_array_order = evecs_array_array


for array_web, array_sat in zip(evecs_array_array, evecs_array_array_sat):
    alpha_array = [] #Ángulo entre lambda_1 de sats y webs.
    beta_array = []
    gamma_array = [] #Ángulo entre lambda_1 de sats y webs.
    delta_array = []
    v1_i_web = evecs_array_array[0][1][0]
    v2_i_web = evecs_array_array[0][1][1]
    v1_i_sat = evecs_array_array_sat[0][1][0]
    v2_i_sat = evecs_array_array_sat[0][1][1]
    alpha_i = Ángulo(v1_i_web,v1_i_sat)
    beta_i = Ángulo(v1_i_web,v2_i_sat)
    del array_web[0]
    del array_sat[0]
    #Se calcula el ángulo entre los autovectores primeros, que corresponden con el menor autovalor:
    for matrix_web, matrix_sat in zip(array_web, array_sat):
        v1_web = matrix_web[0]
        v2_web = matrix_web[1]
        v3_web = matrix_web[2]
        v1_sat = matrix_sat[0]
        v2_sat = matrix_sat[1]
        v3_sat = matrix_sat[2]
        alpha = Ángulo(v1_web,v1_sat)
        beta = Ángulo(v2_web,v2_sat)
        if alpha > math.pi/4:
            alpha = Ángulo(v1_web,v2_sat)
            beta = Ángulo(v2_web,v1_sat)
            alpha_array.append(alpha)
            beta_array.append(beta)
        else:
            alpha_array.append(alpha)
            beta_array.append(beta)
        v1_i_web = v1_web
        v1_i_sat = v1_sat
        v2_i_web = v2_web
        v2_i_sat = v2_sat
            
    fig4 = plt.figure()
    plt.title('Estructura: {} vs {}'.format(webs[m],sats[m]))
    plt.ylabel('Ángulos entre ejes')
    plt.axhline(y=0.785398, xmin=0, xmax=max(atime_array), linewidth = 1)
    plt.axis([min(atime_array), max(atime_array), 0, math.pi/2])
    plt.plot(atime_array, alpha_array, color = "red", linewidth = 1, linestyle = ":")
    plt.plot(atime_array, beta_array, color = "blue", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, alpha_array, color='r', marker='+', s=1e2, linewidths=1, label = r'$\alpha (t)$')
    plt.scatter(atime_array, beta_array, color='b', marker='+', s=1e2, linewidths=1, label = r'$\beta (t)$')
    #plt.plot(atime_array, gamma_array, color = "k", linewidth = 1, linestyle = ":")
    #plt.plot(atime_array, delta_array, color = "y", linewidth = 1, linestyle = ":")
    #plt.scatter(atime_array, gamma_array, color='k', marker='^', s=1e1, linewidths=1, label = r'$\gamma (t)$')
    #plt.scatter(atime_array, delta_array, color='y', marker='^', s=1e1, linewidths=1, label = r'$\delta (t)$')
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    #plt.savefig(images_path + 'Ángulos, atime, eje principal {} vs {}.png'.format(sats[m],webs[m]), bbox_inches='tight', dpi=600)
    plt.savefig(images_path + 'Prueba2 {} vs {}.png'.format(sats[m],webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    m = m + 1

m = 0

for array_web, array_sat in zip(evecs_array_array, evecs_array_array_sat_order):
    alpha_array = [] #Ángulo entre lambda_1 de sats y webs.
    beta_array = []
    gamma_array = []
    v1_i_web = evecs_array_array[0][0][0]
    v2_i_web = evecs_array_array[0][0][1]
    v1_i_sat = evecs_array_array_sat[0][0][0]
    v2_i_sat = evecs_array_array_sat[0][0][1]
    #Se calcula el ángulo entre los autovectores primeros, que corresponden con el menor autovalor:
    for matrix_web, matrix_sat in zip(array_web, array_sat):
        v1_web = matrix_web[0]
        v2_web = matrix_web[1]
        v3_web = matrix_web[2]
        v1_sat = matrix_sat[0]
        v2_sat = matrix_sat[1]
        v3_sat = matrix_sat[2]
        alpha = Ángulo(v1_web,v1_sat)
        beta = Ángulo(v2_web,v2_sat)
        if alpha > math.pi/4:
            alpha = Ángulo(v1_web,v2_sat)
            beta = Ángulo(v2_web,v1_sat)
            alpha_array.append(alpha)
            beta_array.append(beta)
        else:
            alpha_array.append(alpha)
            beta_array.append(beta)
        v1_i_web = v1_web
        v1_i_sat = v1_sat
        v2_i_web = v2_web
        v2_i_sat = v2_sat
            
    fig4 = plt.figure()
    plt.title('Estructura: {} vs {}'.format(webs[m],sats[1 - m]))
    plt.ylabel('Ángulos entre ejes')
    plt.axhline(y=0.785398, xmin=0, xmax=max(atime_array), linewidth = 1)
    plt.axis([min(atime_array), max(atime_array), 0, math.pi/2])
    plt.plot(atime_array, alpha_array, color = "red", linewidth = 1, linestyle = ":")
    plt.plot(atime_array, beta_array, color = "blue", linewidth = 1, linestyle = ":")
    plt.scatter(atime_array, alpha_array, color='r', marker='+', s=1e2, linewidths=1, label = r'$\alpha (t)$')
    plt.scatter(atime_array, beta_array, color='b', marker='+', s=1e2, linewidths=1, label = r'$\beta (t)$')
    #plt.plot(atime_array, gamma_array, color = "k", linewidth = 1, linestyle = ":")
    #plt.plot(atime_array, delta_array, color = "y", linewidth = 1, linestyle = ":")
    #plt.scatter(atime_array, gamma_array, color='k', marker='^', s=1e1, linewidths=1, label = r'$\gamma (t)$')
    #plt.scatter(atime_array, delta_array, color='y', marker='^', s=1e1, linewidths=1, label = r'$\delta (t)$')
    plt.xlabel('t')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
    plt.grid()
    plt.savefig(images_path + 'Ángulos, atime, eje principal {} vs {}.png'.format(sats[1 - m],webs[m]), bbox_inches='tight', dpi=600)
    plt.show()
    
    m = m + 1


fig_COMP = plt.figure()
plt.title('Estructura, comparativa: {} vs {}'.format(webs[0],sats[0]))
plt.ylabel('Cocientes de semiejes, comparativa')
plt.axis([min(atime_array), max(atime_array), 0, 0.6])
plt.plot(atime_array, c_b_array_array[1], color = "red", linewidth = 1, linestyle = ":")
plt.plot(atime_array, c_c_array_array[1], color = "blue", linewidth = 1, linestyle = ":")

plt.plot(atime_array, c_b_array_array_sat[1], color = "red", linewidth = 1, linestyle = "--")
plt.plot(atime_array, c_c_array_array_sat[1], color = "blue", linewidth = 1, linestyle = "--")

line_t = np.linspace(min(atime_array), max(atime_array), 100)

cb_fit_web2 = np.poly1d(np.polyfit(atime_array, c_b_array_array[1], 5))
cc_fit_web2 = np.poly1d(np.polyfit(atime_array, c_c_array_array[1], 5))

cb_fit_sat_no = np.poly1d(np.polyfit(atime_array, c_b_array_array_sat[1], 5))
cc_fit_sat_no = np.poly1d(np.polyfit(atime_array, c_c_array_array_sat[1], 5))


plt.scatter(atime_array, c_b_array_array[1], color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t), web")
plt.scatter(atime_array, c_c_array_array[1], color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t), web")
plt.scatter(atime_array, c_b_array_array_sat[1], color='r', marker='^', s=1e1, linewidths=1, label = "b(t)/a(t), sat")
plt.scatter(atime_array, c_c_array_array_sat[1], color='b', marker='^', s=1e1, linewidths=1, label = "c(t)/a(t), sat")
plt.plot(line_t, cb_fit_web2(line_t),color = 'r', label = 'Ajuste b(t)/a(t), web')
plt.plot(line_t, cc_fit_web2(line_t),color = 'b', label = 'Ajuste c(t)/a(t), web')
plt.plot(line_t, cb_fit_sat_no(line_t),color = 'r', linewidth = 0.5, label = 'Ajuste b(t)/a(t), sat')
plt.plot(line_t, cc_fit_sat_no(line_t),color = 'b', linewidth = 0.5, label = 'Ajuste c(t)/a(t), sat')

plt.grid()
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.savefig(images_path + 'Semiejes, atime, cocientes, comparativa {} vs {}.png'.format(webs[0],sats[0]), bbox_inches='tight', dpi=600)
plt.show()


fig_COMP = plt.figure()
plt.title('Estructura, comparativa: {} vs {}'.format(webs[0],sats[1]))
plt.ylabel('Cocientes de semiejes, comparativa')
plt.axis([min(atime_array), max(atime_array), 0, 0.6])
plt.plot(atime_array, c_b_array_array[1], color = "red", linewidth = 1, linestyle = ":")
plt.plot(atime_array, c_c_array_array[1], color = "blue", linewidth = 1, linestyle = ":")

plt.plot(atime_array, c_b_array_array_sat[3], color = "red", linewidth = 1, linestyle = "--")
plt.plot(atime_array, c_c_array_array_sat[3], color = "blue", linewidth = 1, linestyle = "--")

line_t = np.linspace(min(atime_array), max(atime_array), 100)

cb_fit_web2 = np.poly1d(np.polyfit(atime_array, c_b_array_array[1], 5))
cc_fit_web2 = np.poly1d(np.polyfit(atime_array, c_c_array_array[1], 5))

cb_fit_sat = np.poly1d(np.polyfit(atime_array, c_b_array_array_sat[3], 5))
cc_fit_sat = np.poly1d(np.polyfit(atime_array, c_c_array_array_sat[3], 5))


plt.scatter(atime_array, c_b_array_array[1], color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t), web")
plt.scatter(atime_array, c_c_array_array[1], color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t), web")
plt.scatter(atime_array, c_b_array_array_sat[3], color='r', marker='^', s=1e1, linewidths=1, label = "b(t)/a(t), sat")
plt.scatter(atime_array, c_c_array_array_sat[3], color='b', marker='^', s=1e1, linewidths=1, label = "c(t)/a(t), sat")
plt.plot(line_t, cb_fit_web2(line_t),color = 'r', label = 'Ajuste b(t)/a(t), web')
plt.plot(line_t, cc_fit_web2(line_t),color = 'b', label = 'Ajuste c(t)/a(t), web')
plt.plot(line_t, cb_fit_sat(line_t),color = 'r', linewidth = 0.5, label = 'Ajuste b(t)/a(t), sat')
plt.plot(line_t, cc_fit_sat(line_t),color = 'b', linewidth = 0.5, label = 'Ajuste c(t)/a(t), sat')

plt.grid()
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.savefig(images_path + 'Semiejes, atime, cocientes, comparativa {} vs {}.png'.format(webs[0],sats[1]), bbox_inches='tight', dpi=600)
plt.show()


fig_COMP = plt.figure()
plt.title('Estructura, comparativa: {} vs {}'.format(webs[1],sats[0]))
plt.ylabel('Cocientes de semiejes, comparativa')
plt.axis([min(atime_array), max(atime_array), 0, 0.6])
plt.plot(atime_array, c_b_array_array[3], color = "red", linewidth = 1, linestyle = ":")
plt.plot(atime_array, c_c_array_array[3], color = "blue", linewidth = 1, linestyle = ":")

plt.plot(atime_array, c_b_array_array_sat[1], color = "red", linewidth = 1, linestyle = "--")
plt.plot(atime_array, c_c_array_array_sat[1], color = "blue", linewidth = 1, linestyle = "--")

line_t = np.linspace(min(atime_array), max(atime_array), 100)

cb_fit_web3 = np.poly1d(np.polyfit(atime_array, c_b_array_array[3], 5))
cc_fit_web3 = np.poly1d(np.polyfit(atime_array, c_c_array_array[3], 5))

cb_fit_sat_no = np.poly1d(np.polyfit(atime_array, c_b_array_array_sat[1], 5))
cc_fit_sat_no = np.poly1d(np.polyfit(atime_array, c_c_array_array_sat[1], 5))


plt.scatter(atime_array, c_b_array_array[3], color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t), web")
plt.scatter(atime_array, c_c_array_array[3], color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t), web")
plt.scatter(atime_array, c_b_array_array_sat[1], color='r', marker='^', s=1e1, linewidths=1, label = "b(t)/a(t), sat")
plt.scatter(atime_array, c_c_array_array_sat[1], color='b', marker='^', s=1e1, linewidths=1, label = "c(t)/a(t), sat")
plt.plot(line_t, cb_fit_web3(line_t),color = 'r', label = 'Ajuste b(t)/a(t), web')
plt.plot(line_t, cc_fit_web3(line_t),color = 'b', label = 'Ajuste c(t)/a(t), web')
plt.plot(line_t, cb_fit_sat_no(line_t),color = 'r', linewidth = 0.5, label = 'Ajuste b(t)/a(t), sat')
plt.plot(line_t, cc_fit_sat_no(line_t),color = 'b', linewidth = 0.5, label = 'Ajuste c(t)/a(t), sat')

plt.grid()
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.savefig(images_path + 'Semiejes, atime, cocientes, comparativa {} vs {}.png'.format(webs[1],sats[0]), bbox_inches='tight', dpi=600)
plt.show()


fig_COMP = plt.figure()
plt.title('Estructura, comparativa: {} vs {}'.format(webs[1],sats[1]))
plt.ylabel('Cocientes de semiejes, comparativa')
plt.axis([min(atime_array), max(atime_array), 0, 0.6])
plt.plot(atime_array, c_b_array_array[3], color = "red", linewidth = 1, linestyle = ":")
plt.plot(atime_array, c_c_array_array[3], color = "blue", linewidth = 1, linestyle = ":")

plt.plot(atime_array, c_b_array_array_sat[3], color = "red", linewidth = 1, linestyle = "--")
plt.plot(atime_array, c_c_array_array_sat[3], color = "blue", linewidth = 1, linestyle = "--")

line_t = np.linspace(min(atime_array), max(atime_array), 100)

cb_fit_web3 = np.poly1d(np.polyfit(atime_array, c_b_array_array[3], 5))
cc_fit_web3 = np.poly1d(np.polyfit(atime_array, c_c_array_array[3], 5))

cb_fit_sat = np.poly1d(np.polyfit(atime_array, c_b_array_array_sat[3], 5))
cc_fit_sat = np.poly1d(np.polyfit(atime_array, c_c_array_array_sat[3], 5))


plt.scatter(atime_array, c_b_array_array[3], color='r', marker='+', s=1e2, linewidths=1, label = "b(t)/a(t), web")
plt.scatter(atime_array, c_c_array_array[3], color='b', marker='+', s=1e2, linewidths=1, label = "c(t)/a(t), web")
plt.scatter(atime_array, c_b_array_array_sat[3], color='r', marker='^', s=1e1, linewidths=1, label = "b(t)/a(t), sat")
plt.scatter(atime_array, c_c_array_array_sat[3], color='b', marker='^', s=1e1, linewidths=1, label = "c(t)/a(t), sat")
plt.plot(line_t, cb_fit_web3(line_t),color = 'r', label = 'Ajuste b(t)/a(t), web')
plt.plot(line_t, cc_fit_web3(line_t),color = 'b', label = 'Ajuste c(t)/a(t), web')
plt.plot(line_t, cb_fit_sat(line_t),color = 'r', linewidth = 0.5, label = 'Ajuste b(t)/a(t), sat')
plt.plot(line_t, cc_fit_sat(line_t),color = 'b', linewidth = 0.5, label = 'Ajuste c(t)/a(t), sat')

plt.grid()
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.savefig(images_path + 'Semiejes, atime, cocientes, comparativa {} vs {}.png'.format(webs[1],sats[1]), bbox_inches='tight', dpi=600)
plt.show()


figp1_COMP = plt.figure()
plt.title('Triaxialidad, comparativa: {} vs {}'.format(webs[0],sats[0]))
plt.ylabel('T(t)')
plt.axis([min(atime_array), max(atime_array), min(T_sat_array_array[1])-min(T_sat_array_array[1])*0.1, max(T_web_array_array[1])+max(T_web_array_array[1])*0.1])
plt.plot(atime_array, T_web_array_array[1], color = "red", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, T_web_array_array[1], color='r', marker='+', s=1e2, linewidths=1, label = "T(t), web")
plt.plot(atime_array, T_sat_array_array[1], color = "red", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, T_sat_array_array[1], color='r', marker='^', s=5e1, linewidths=1, label = "T(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Triaxialidad, comparativa {} vs {}.png'.format(webs[0],sats[0]), bbox_inches='tight', dpi=600)
plt.show()

figp2_COMP = plt.figure()
plt.title('Elipticidad, comparativa: {} vs {}'.format(webs[0],sats[0]))
plt.ylabel('e(t)')
plt.axis([min(atime_array), max(atime_array), min(e_sat_array_array[1])-min(e_sat_array_array[1])*0.1, max(e_web_array_array[1])+max(e_web_array_array[1])*0.1])
plt.plot(atime_array, e_web_array_array[1], color = "b", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, e_web_array_array[1], color='b', marker='+', s=1e2, linewidths=1, label = "e(t), web")
plt.plot(atime_array, e_sat_array_array[1], color = "b", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, e_sat_array_array[1], color='b', marker='^', s=5e1, linewidths=1, label = "e(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Elipticidad, comparativa {} vs {}.png'.format(webs[0],sats[0]), bbox_inches='tight', dpi=600)
plt.show()

figp3_COMP = plt.figure()
plt.title('Prolaticidad, comparativa: {} vs {}'.format(webs[0],sats[0]))
plt.ylabel('p(t)')
plt.axis([min(atime_array), max(atime_array), min(p_sat_array_array[1])-min(p_sat_array_array[1])*0.1, max(p_web_array_array[1])+max(p_web_array_array[1])*0.1])
plt.plot(atime_array, p_web_array_array[1], color = "g", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, p_web_array_array[1], color='g', marker='+', s=1e2, linewidths=1, label = "p(t), web")
plt.plot(atime_array, p_sat_array_array[1], color = "g", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, p_sat_array_array[1], color='g', marker='^', s=5e1, linewidths=1, label = "p(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Prolaticidad, comparativa {} vs {}.png'.format(webs[0],sats[0]), bbox_inches='tight', dpi=600)
plt.show()


figp1_COMP = plt.figure()
plt.title('Triaxialidad, comparativa: {} vs {}'.format(webs[0],sats[1]))
plt.ylabel('T(t)')
plt.axis([min(atime_array), max(atime_array), min(T_sat_array_array[1])-min(T_sat_array_array[1])*0.1, max(T_web_array_array[1])+max(T_web_array_array[1])*0.1])
plt.plot(atime_array, T_web_array_array[1], color = "red", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, T_web_array_array[1], color='r', marker='+', s=1e2, linewidths=1, label = "T(t), web")
plt.plot(atime_array, T_sat_array_array[3], color = "red", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, T_sat_array_array[3], color='r', marker='^', s=5e1, linewidths=1, label = "T(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Triaxialidad, comparativa {} vs {}.png'.format(webs[0],sats[1]), bbox_inches='tight', dpi=600)
plt.show()

figp2_COMP = plt.figure()
plt.title('Elipticidad, comparativa: {} vs {}'.format(webs[0],sats[1]))
plt.ylabel('e(t)')
plt.axis([min(atime_array), max(atime_array), min(e_sat_array_array[1])-min(e_sat_array_array[1])*0.1, max(e_web_array_array[1])+max(e_web_array_array[1])*0.1])
plt.plot(atime_array, e_web_array_array[1], color = "b", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, e_web_array_array[1], color='b', marker='+', s=1e2, linewidths=1, label = "e(t), web")
plt.plot(atime_array, e_sat_array_array[3], color = "b", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, e_sat_array_array[3], color='b', marker='^', s=5e1, linewidths=1, label = "e(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Elipticidad, comparativa {} vs {}.png'.format(webs[0],sats[1]), bbox_inches='tight', dpi=600)
plt.show()

figp3_COMP = plt.figure()
plt.title('Prolaticidad, comparativa: {} vs {}'.format(webs[0],sats[1]))
plt.ylabel('p(t)')
plt.axis([min(atime_array), max(atime_array), min(p_sat_array_array[1])-min(p_sat_array_array[1])*0.1, max(p_web_array_array[1])+max(p_web_array_array[1])*0.1])
plt.plot(atime_array, p_web_array_array[1], color = "g", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, p_web_array_array[1], color='g', marker='+', s=1e2, linewidths=1, label = "p(t), web")
plt.plot(atime_array, p_sat_array_array[3], color = "g", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, p_sat_array_array[3], color='g', marker='^', s=5e1, linewidths=1, label = "p(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Prolaticidad, comparativa {} vs {}.png'.format(webs[0],sats[1]), bbox_inches='tight', dpi=600)
plt.show()


figp1_COMP = plt.figure()
plt.title('Triaxialidad, comparativa: {} vs {}'.format(webs[1],sats[0]))
plt.ylabel('T(t)')
plt.axis([min(atime_array), max(atime_array), min(T_sat_array_array[3])-min(T_sat_array_array[3])*0.1, max(T_web_array_array[3])+max(T_web_array_array[3])*0.1])
plt.plot(atime_array, T_web_array_array[3], color = "red", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, T_web_array_array[3], color='r', marker='+', s=1e2, linewidths=1, label = "T(t), web")
plt.plot(atime_array, T_sat_array_array[1], color = "red", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, T_sat_array_array[1], color='r', marker='^', s=5e1, linewidths=1, label = "T(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Triaxialidad, comparativa {} vs {}.png'.format(webs[1],sats[0]), bbox_inches='tight', dpi=600)
plt.show()

figp2_COMP = plt.figure()
plt.title('Elipticidad, comparativa: {} vs {}'.format(webs[1],sats[0]))
plt.ylabel('e(t)')
plt.axis([min(atime_array), max(atime_array), min(e_sat_array_array[3])-min(e_sat_array_array[3])*0.1, max(e_web_array_array[3])+max(e_web_array_array[3])*0.1])
plt.plot(atime_array, e_web_array_array[3], color = "b", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, e_web_array_array[3], color='b', marker='+', s=1e2, linewidths=1, label = "e(t), web")
plt.plot(atime_array, e_sat_array_array[1], color = "b", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, e_sat_array_array[1], color='b', marker='^', s=5e1, linewidths=1, label = "e(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Elipticidad, comparativa {} vs {}.png'.format(webs[1],sats[0]), bbox_inches='tight', dpi=600)
plt.show()

figp3_COMP = plt.figure()
plt.title('Prolaticidad, comparativa: {} vs {}'.format(webs[1],sats[0]))
plt.ylabel('p(t)')
plt.axis([min(atime_array), max(atime_array), min(p_sat_array_array[3])-min(p_sat_array_array[3])*0.1, max(p_web_array_array[3])+max(p_web_array_array[3])*0.1])
plt.plot(atime_array, p_web_array_array[3], color = "g", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, p_web_array_array[3], color='g', marker='+', s=1e2, linewidths=1, label = "p(t), web")
plt.plot(atime_array, p_sat_array_array[1], color = "g", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, p_sat_array_array[1], color='g', marker='^', s=5e1, linewidths=1, label = "p(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Prolaticidad, comparativa {} vs {}.png'.format(webs[1],sats[0]), bbox_inches='tight', dpi=600)
plt.show()


figp1_COMP = plt.figure()
plt.title('Triaxialidad, comparativa: {} vs {}'.format(webs[1],sats[1]))
plt.ylabel('T(t)')
plt.axis([min(atime_array), max(atime_array), min(T_sat_array_array[3])-min(T_sat_array_array[3])*0.1, max(T_web_array_array[3])+max(T_web_array_array[3])*0.1])
plt.plot(atime_array, T_web_array_array[3], color = "red", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, T_web_array_array[3], color='r', marker='+', s=1e2, linewidths=1, label = "T(t), web")
plt.plot(atime_array, T_sat_array_array[3], color = "red", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, T_sat_array_array[3], color='r', marker='^', s=5e1, linewidths=1, label = "T(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Triaxialidad, comparativa {} vs {}.png'.format(webs[1],sats[1]), bbox_inches='tight', dpi=600)
plt.show()

figp2_COMP = plt.figure()
plt.title('Elipticidad, comparativa: {} vs {}'.format(webs[1],sats[1]))
plt.ylabel('e(t)')
plt.axis([min(atime_array), max(atime_array), min(e_sat_array_array[3])-min(e_sat_array_array[3])*0.1, max(e_web_array_array[3])+max(e_web_array_array[3])*0.1])
plt.plot(atime_array, e_web_array_array[3], color = "b", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, e_web_array_array[3], color='b', marker='+', s=1e2, linewidths=1, label = "e(t), web")
plt.plot(atime_array, e_sat_array_array[3], color = "b", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, e_sat_array_array[3], color='b', marker='^', s=5e1, linewidths=1, label = "e(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Elipticidad, comparativa {} vs {}.png'.format(webs[1],sats[1]), bbox_inches='tight', dpi=600)
plt.show()

figp3_COMP = plt.figure()
plt.title('Prolaticidad, comparativa: {} vs {}'.format(webs[1],sats[1]))
plt.ylabel('p(t)')
plt.axis([min(atime_array), max(atime_array), min(p_sat_array_array[3])-min(p_sat_array_array[3])*0.1, max(p_web_array_array[3])+max(p_web_array_array[3])*0.1])
plt.plot(atime_array, p_web_array_array[3], color = "g", linewidth = 1, linestyle = ":")
plt.scatter(atime_array, p_web_array_array[3], color='g', marker='+', s=1e2, linewidths=1, label = "p(t), web")
plt.plot(atime_array, p_sat_array_array[3], color = "g", linewidth = 1, linestyle = "--")
plt.scatter(atime_array, p_sat_array_array[3], color='g', marker='^', s=5e1, linewidths=1, label = "p(t), sat")
plt.xlabel('t')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=5)
plt.grid()
plt.savefig(images_path + 'Prolaticidad, comparativa {} vs {}.png'.format(webs[1],sats[1]), bbox_inches='tight', dpi=600)
plt.show()

end = time.time()
print('___________________________________________________________________________________________________________')
print("Tiempo de ejecución (s):", end-start)
print("En minutos (min):", (end-start)/60)
print('___________________________________________________________________________________________________________')



