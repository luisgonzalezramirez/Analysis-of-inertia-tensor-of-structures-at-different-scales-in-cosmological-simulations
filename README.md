# Analysis-of-inertia-tensor-of-structures-at-different-scales-in-cosmological-simulations

## Astrophysics Department
###### Luis González Ramírez

_________________________________________________

Physics final degree project, Universidad Complutense de Madrid.

The codes used for the analysis of the inertial tensor of cosmological structures of an N-body simulation are presented here. The main idea of sharing these codes is to allow the reproduction of the results (inertial tensor, eigenvectors, evolution animations, centers of mass, etc.) for simulations of this type. However, it is not possible to share the data used due to copyright. Even so, these codes can be useful as a guide.


The framework of the computational method that make up these codes has been designed in such a way that enable to run each code individually, except those codes designated with (sat), that need to run the code Satellite_reading and Satellite_concatenation previously, its structured as follows:

1º/ Satellites classification and particles selection:
  - Satellite_reading, reads the two list of satellite spheres, and select particles of each one.
  - Satellite_concatenation, gather all satellite particles into two .txt files.

2º/ Parameters reading and overview of simulation (types of particles, structures and lagrangian volumes):
  - Parameters_reading.py, gets simulation initial parameters and conversion units.
  - Parameters_reading_and_simulation_box_size.py, shows simulation box size and cosmological parameters
  - Parameters_reading_and_simulation_frames, graphs every simulation frame (plotting distinctively
 all particles types) along the conversion and parameters units.
 
 3º/ Calculations and results:
  - Center_of_mass.
  - Convergence_of_web_structures, studies the convergence of web structures center of mass.
  - Inertia_tensor, all calculations of axis, eigenvectors, tensor components, morphologic parameters and angles between axis.
  - Eigenvectors_axis_satellites_and_webs, plots eigenvectors of each pair of web and satellite group.

Making GIF simulations (for Parameters_reading_and_simulation_frames and Eigenvectors_axis_satellites_and_webs) are done with two codes: SimulationGIF and EigenvectorsGIF
