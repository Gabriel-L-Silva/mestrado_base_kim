import cv2
import numpy as np
import glob
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
img_array = []
for filename in glob.glob('./../sim/*'):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)


out = cv2.VideoWriter('./../sim_video.avi',cv2.VideoWriter_fourcc(*'DIVX'), 15, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()