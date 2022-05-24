import glob

from PIL import Image
 
frames = []
imgs = glob.glob("*.png")
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
 

frames[0].save('png_to_gif.gif', format='GIF', append_images=frames[1:], save_all=True, duration=300, loop=0)