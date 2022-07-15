import cv2
import numpy as np
import glob
import os


def number(filename):
    word, rest = filename.split('_')
    number, ext = rest.split('.')
    return (word, int(number))

dir_path = os.path.dirname(os.path.realpath(__file__))
img_array = []
images = sorted(glob.glob('./../sim/*'), key=number)
for filename in images:
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)


out = cv2.VideoWriter('./../sim_video.avi',cv2.VideoWriter_fourcc(*'DIVX'), 60, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()