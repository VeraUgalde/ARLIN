# -*- coding: utf-8 -*-
"""Neuro_Imaging_pt1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1vx8ZG1UvZLwICwjdIV5EP6mxn29KpHNV
"""

import os
import glob
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize, disk, binary_dilation, binary_erosion
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage import measure
from datetime import date
from fil_finder import FilFinder2D
import astropy.units as u

def main():
    #Getting path to Images
    print("\n Full path to folder of images: ")
    directory = input()
    
    while not os.path.isdir(directory):
        print("\n Incorrect path, please try again:")
        print("\n Full path to folder of images: ")
        directory = input()
    
    os.chdir(directory)
    
    #Making new subfolder to store skeletons and prints
    print_dir = "SomaAndNucleusPrints"
    if not os.path.exists(print_dir):
        os.makedirs(print_dir)  
    
    #Getting name of DAPI channel
    dapi_channel = " " 
    print("\n Name of DAP channel: ")
    print("\n ex: CY7, c7, DAPI, etc")
    dapi_channel = input()
    while not glob.glob("*XY*" + dapi_channel + ".tif"):
        print("\n Incorrect channel name. Try again: ")
        print("\n Name of DAPI channel: ")
        print("\n ex: CY7, c7, DAPI, etc")
        dapi_channel = input()
    
    
    #Getting mRNA channels
    mRNA_channel = " "
    mRNA_channels = []
    while mRNA_channel != "0":
        print("\n Name of mRNA channel (if done press 0): ")
        print("\n ex: CY5, c5, CY3, c3 etc")
        mRNA_channel = input()
        while not glob.glob("*XY*" + mRNA_channel + ".tif") and mRNA_channel != "0": 
            print("\n Incorrect channel name. Try again: ")
            print("\n Name of mRNA channel (if done press 0): ")
            print("\n ex: CY5, c5, CY3, c3 etc")
            mRNA_channel = input()
        if mRNA_channel != "0":
            mRNA_channels.append(mRNA_channel)
    
    #Getting MAP2 channel
    map2_channel = " " 
    print("\n Name of MAP2 channel: ")
    print("\n ex: CY7, c7, DAPI, etc")
    map2_channel = input()
    while not glob.glob("*XY*" + map2_channel + ".tif"):
        print("\n Incorrect channel name. Try again: ")
        print("\n Name of MAP2 channel: ")
        print("\n ex: CY7, c7, DAPI, etc")
        map2_channel = input()
    
    
    
    max_files = glob.glob("MAX*" + map2_channel + ".gif")
    max_files = ["MAX_Exp36_DIV14_MG6R_HSPA8MAP2PSD95DAPI_024_XY1_CY7.gif"]
    for max_file in max_files:
        print(os.path.basename(max_file)[:-4])
        color_file = glob.glob(os.path.basename(max_file)[:-4] + "*color.gif")
        
        if len(color_file) > 0:
            color_file = color_file[0]
            channel_index = max_file.index(map2_channel)
            dapi_file = glob.glob(max_file[:channel_index] + dapi_channel + "*.gif")[0]
        
            soma_prints, nuc_prints = process_img(dapi_file, color_file)
    
            #Save skeletons in "/processed/FOV#" folder
            for i in range(len(soma_prints)):
                Image.fromarray(soma_prints[i]).save(print_dir +  "/" + os.path.basename(max_file)[:-4] + "_soma_print_" + str(i+1) + ".gif")
                Image.fromarray(nuc_prints[i]).save(print_dir +  "/" + os.path.basename(max_file)[:-4] + "_nuc_print_" + str(i+1) + ".gif")
              
            for mRNA_channel in mRNA_channels:
                end_index = max_file.index(map2_channel)
                print(os.path.basename(max_file)[4:end_index] + mRNA_channel + ".tif")
                raw_file = glob.glob(os.path.basename(max_file)[4:end_index] + mRNA_channel + ".tif")[0]
                dapi_file = glob.glob(os.path.basename(max_file)[4:end_index] + dapi_channel + ".tif")[0]
                  
                #Making new subfolder to store outlines each channel
                channel_dir = "Outlines_" + mRNA_channel
                if not os.path.exists(channel_dir):
                  os.makedirs(channel_dir) 
                
                
                
                getOutlinesTxts(raw_file, dapi_file, channel_dir, soma_prints, nuc_prints)
    
    return

def process_img(dapi_file, color_file):

  #reading in colored image and regular image
  color_img = plt.imread(color_file)[:,:,0:3] #should be colored map2
  dapi_img = plt.imread(dapi_file) #should be dapi
  
  #extracting the number of colors used to mark dendrites in color_img
  num_colors = 0
  i = color_file.index("_color.gif")
  if color_file[i-2] == "_":
    num_colors = int(color_file[i-1])
  else:
    num_colors = int(color_file[i-2:i])

  #Pretty colors yay!!
  red  = [255,0,0]
  green = [0,255,0]
  blue   = [0,0,255]
  orange = [255, 147, 0]
  yellow = [255, 255, 0]
  purple = [153, 41, 189]
  teal  = [83, 219, 196]
  mint = [204, 252, 216]
  maroon = [159, 17, 0]
  salmon = [255, 128, 110]
  colors = [red, green, blue, orange, yellow, purple, teal, mint, maroon, salmon]


  soma_prints = [None]*num_colors
  nuc_prints = [None]*num_colors
  
  #extracting individual dendrites from image
  j = 0
  for color in colors:
    if j < num_colors:
      mask = np.where(np.all(color_img == color, axis=-1),1, 0)
      nucleus = np.multiply(dapi_img, mask)
      mask = np.where(mask==1, 255, 0)
      soma_prints[j] = mask.astype(np.uint8)
      nuc_prints[j] = getNucleus(nucleus, dapi_img)

      j+= 1
  
  return soma_prints, nuc_prints

def getNucleus(nucleus, img):
  
  #initial thresholding using whole image
  thresh = threshold_otsu(img) * 1.2
  
  nuc_print = nucleus > thresh
  return nuc_print


def getOutlinesTxts(channel_img, dapi_img, channel_dir, soma_prints, nuc_prints):

  if len(soma_prints) == 1:
    soma_contours = measure.find_contours(soma_prints[0], 0)
    soma_contours = np.array(soma_contours)

    nuc_contours = measure.find_contours(nuc_prints[0], 0)
    nuc_contours = np.array(nuc_contours)
    if np.shape(soma_contours)[0] == 0 or np.shape(nuc_contours)[0] == 0:
        return
      
  outline_txt = open(channel_dir + "/" + channel_img[:-4] + "__outline.txt", "a")

  outline_txt.write("FISH-QUANT\n")
  outline_txt.write("File-version   3D_v1\n")
  outline_txt.write("RESULTS OF SPOT DETECTION PERFORMED ON " + date.today().strftime("%d-%b-%Y") + "\n")
  outline_txt.write("COMMENT	Automated outline definition (batch or quick-save)\n")
  outline_txt.write("IMG_Raw\t" + channel_img + "\n")
  outline_txt.write("IMG_Filtered\n")
  outline_txt.write("IMG_DAPI\t" + dapi_img + "\n")
  outline_txt.write("IMG_TS_label\n")
  outline_txt.write("FILE_settings\n")
  outline_txt.write("PARAMETERS\n")
  outline_txt.write("Pix-XY	Pix-Z	RI	Ex	Em	NA	Type\n")
  outline_txt.write("107.5	300	1.35	547	583	1.4	widefield\n")

  for i in range(len(soma_prints)):

    #getting counterclockwise list of pixels on outlines
    soma_contours = measure.find_contours(soma_prints[i], 0)
    soma_contours = np.array(soma_contours)

    nuc_contours = measure.find_contours(nuc_prints[i], 0)
    nuc_contours = np.array(nuc_contours)

    if np.shape(nuc_contours)[0] >= 2 or np.shape(soma_contours)[0] >= 2:
        print("This soma is sketchy")
        soma_contours = np.array([find_max_list(soma_contours)])
        nuc_contours = np.array([find_max_list(nuc_contours)])
    elif np.shape(soma_contours)[0] == 0:
        print("Colours are wrong in this image!")
        break
    

    soma_contours = soma_contours.reshape(soma_contours.shape[1:])
    nuc_contours = nuc_contours.reshape(nuc_contours.shape[1:])

    outline_txt.write("CELL_START\tCell_" + str(i+1) + "\n")

    outline_txt.write("X_POS")
    for i in range(len(soma_contours)):
     if i % 10 == 0:
        outline_txt.write("\t" + str(int(soma_contours[i][1])))
    outline_txt.write("\t")

    outline_txt.write("\nY_POS")
    for i in range(len(soma_contours)):
      if i % 10 == 0:
        outline_txt.write("\t" + str(int(soma_contours[i][0])))
    outline_txt.write("\t")

    outline_txt.write("\nZ_Pos\t\n")
    outline_txt.write("CELL_END\n")
    outline_txt.write("Nucleus_START\tNucleus_manual\n")

    outline_txt.write("X_POS")
    for i in range(len(nuc_contours)):
     if i % 10 == 0:
        outline_txt.write("\t" + str(int(nuc_contours[i][1])))
    outline_txt.write("\t")

    outline_txt.write("\nY_POS")
    for i in range(len(nuc_contours)):
      if i % 10 == 0:
        outline_txt.write("\t" + str(int(nuc_contours[i][0])))
    outline_txt.write("\t")

    outline_txt.write("\nZ_Pos\t\n")
    outline_txt.write("Nucleus_END\n")

    print("made it through")
  outline_txt.close()

  return


def find_max_list(list):
    list_len = [len(i) for i in list]
    return list[np.argmax(np.array(list_len))]

main()
