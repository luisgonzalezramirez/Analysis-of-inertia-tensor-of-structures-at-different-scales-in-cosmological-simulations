# Analysis of inertia tensor of structures at different scales in cosmological simulations

## Astrophysics Department
###### Luis González Ramírez

_________________________________________________

Physics final degree project, Universidad Complutense de Madrid.

The codes used for the analysis of the inertial tensor of cosmological structures of an N-body simulation are presented here. The main idea of sharing these codes is to allow the reproduction of the results (inertial tensor, eigenvectors, evolution animations, centers of mass, etc.) for simulations of this type. However, it is not possible to share the data used due to copyright. Even so, these codes can be useful as a guide.


The framework of the computational method that make up these codes has been designed in such a way that enable to run each code individually, except those codes designated with (sat), that need to run the code Satellite_reading.py and Satellite_concatenation.py previously, its structured as follows:

1º/ Satellites classification and particles selection:
  - Satellite_reading.py, reads the two list of satellite spheres, and select particles of each one.
  - Satellite_concatenation.py, gather all satellite particles into two .txt files.

2º/ Parameters reading and overview of simulation (types of particles, structures and lagrangian volumes):
  - Parameters_reading.py, gets simulation initial parameters and conversion units.
  - Parameters_reading_and_simulation_box_size.py, shows simulation box size and cosmological parameters
  - Parameters_reading_and_simulation_frames.py, graphs every simulation frame along the redshift values (plotting distinctively
 all particles types) beside the conversion and parameters units.
 
 3º/ Calculations and results:
  - Center_of_mass.py (sat).
  - Convergence_of_web_structures.py, studies the convergence of web structures center of mass.
  - Inertia_tensor.py (sat), all calculations of axis, eigenvectors, tensor components, morphologic parameters and angles between axis.
  - Eigenvectors_axis_satellites_and_webs.py (sat), plots eigenvectors of each pair of web and satellite group.

Making GIF simulations (for Parameters_reading_and_simulation_frames.py and Eigenvectors_axis_satellites_and_webs.py) are done with two codes: SimulationGIF.py and EigenvectorsGIF.py.

-------------------------------------------------
Every code still has annotations in Spanish. In further updates they will be changed.
-------------------------------------------------

_________________________________________________
_________________________________________________

# Análisis del tensor de inercia de estructuras a diferentes escalas en simulaciones cosmológicas

## Departamento de Astrofísica
###### Luis González Ramírez

_________________________________________________

Trabajo final de carrera de Física, Universidad Complutense de Madrid.

Aquí se presentan los códigos utilizados para el análisis del tensor de inercia de estructuras cosmológicas de una simulación de N-cuerpos. La idea principal de compartir estos códigos es permitir la reproducción de los resultados (tensor de inercia, vectores propios, animaciones de evolución, centros de masa, etc.) para simulaciones de este tipo. Sin embargo, no es posible compartir los datos utilizados debido a los derechos de autor. Aun así, estos códigos pueden ser útiles como guía.


El marco del método computacional que componen estos códigos ha sido diseñado de tal manera que permite ejecutar cada código individualmente, excepto aquellos códigos designados con (sat), que necesitan ejecutar el código Satellite_reading.py y Satellite_concatenation.py previamente, está estructurado de la siguiente manera:

1º/ Clasificación de satélites y selección de partículas:
  - Satellite_reading.py, lee la lista de dos esferas satelitales y selecciona partículas de cada una.
  - Satellite_concatenation.py, reúne todas las partículas de satélite en dos archivos .txt.

2º/ Lectura de parámetros y resumen de simulación (tipos de partículas, estructuras y volúmenes lagrangianos):
  - Parameters_reading.py, obtiene los parámetros iniciales de la simulación y las unidades de conversión.
  - Parameters_reading_and_simulation_box_size.py, muestra el tamaño del cuadro de simulación y los parámetros cosmológicos
  - Parameters_reading_and_simulation_frames.py, grafica cada cuadro de simulación a lo largo de los valores de corrimiento al rojo (trazando distintivamente
 todos los tipos de partículas) a lo largo de las unidades de conversión y parámetro.
 
 3º/ Cálculos y resultados:
  - Center_of_mass.py (sat).
  - Convergence_of_web_structures.py, estudia la convergencia del centro de masa de las estructuras web.
  - Inertia_tensor.py (sat), todos los cálculos de eje, vectores propios, componentes del tensor, parámetros morfológicos y ángulos entre ejes.
  - Eigenvectors_axis_satellites_and_webs.py (sat), traza los vectores propios de cada par de web y grupo de satélites.

Las simulaciones GIF (para Parameters_reading_and_simulation_frames.py y Eigenvectors_axis_satellites_and_webs.py) se realizan con dos códigos: SimulationGIF.py y EigenvectorsGIF.py.
 
