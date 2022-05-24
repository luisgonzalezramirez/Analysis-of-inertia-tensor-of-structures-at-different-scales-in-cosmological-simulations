
#Función para la concatenación de las partículas de los satélites:
#___________________________________________________________#


from os import listdir

print('___________________________________________________________________________________________________________')
target = input('Introduce dirección local del código:\n')

print('___________________________________________________________________________________________________________')
satellites_particles_path = input('Introduce dirección local de entrada de las listas de partículas de SP:\n')

print('___________________________________________________________________________________________________________')
no_satellites_particles_path = input('Introduce dirección local de entrada de las listas de partículas de SnP:\n')

sat_txt = listdir(satellites_particles_path)
sat_no_txt = listdir(no_satellites_particles_path)

def Concatenación_Satélites(sat_txt,sat_no_txt):
    
    '''
    Función:
    ________________________________________________________________________________________
    Concatena todos los archivos .txt de los satélites .
    
    '''
    
    with open('Partículas No Persistentes.txt','w') as outfile:
        for filename in sat_no_txt:
            with open(filename) as infile:
                contents = infile.read()
                outfile.write(contents)
            
    with open('Partículas Persistentes.txt','w') as outfile:
        for filename in sat_txt:
            with open(filename) as infile:
                contents = infile.read()
                outfile.write(contents)

Concatenación_Satélites(sat_txt,sat_no_txt)