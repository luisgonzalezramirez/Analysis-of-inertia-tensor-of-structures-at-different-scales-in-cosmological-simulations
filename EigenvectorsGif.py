
#Código de creación de un GIF de los ejes:
#_______________________________________________________________#

#Este código se implementa en cada carpeta de estructura web para obtener un GIF.

#Importar módulos de libreria de python:

from PIL import Image
import glob
import time

start = time.time()

frames = []
imgs = glob.glob("*.png")
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)

frames[0].save('Simulación Red C.gif', format='GIF',append_images=frames[1:],save_all=True,duration=300, loop=0)


end = time.time()
print('___________________________________________________________________________________________________________')
print("Tiempo de ejecución (s):", end-start)
print("En minutos (min):", (end-start)/60)
print('___________________________________________________________________________________________________________')

