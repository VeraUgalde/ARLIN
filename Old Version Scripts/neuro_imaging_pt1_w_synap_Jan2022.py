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
from matplotlib.pyplot import imread
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
    skel_dir = "SkeletonsAndPrints"
    if not os.path.exists(skel_dir):
        os.makedirs(skel_dir)  
    
    #Getting name of DAPI channel
    dapi_channel = " " 
    print("\n Name of DAPI channel: ")
    print("\n ex: CY7, c7, DAPI, etc")
    dapi_channel = input()
    while not glob.glob("*XY*" + dapi_channel + ".tif"):
        print("\n Incorrect channel name. Try again: ")
        print("\n Name of DAPI channel: ")
        print("\n ex: CY7, c7, DAPI, etc")
        dapi_channel = input()

    #Getting name of syanpse channel
    synapse_channel = " " 
    print("\n Name of synapse channel: ")
    print("\n ex: CY7, c7, DAPI, etc")
    synapse_channel = input()
    while not glob.glob("*XY*" + synapse_channel + ".tif"):
        print("\n Incorrect channel name. Try again: ")
        print("\n Name of synapse channel: ")
        print("\n ex: CY7, c7, DAPI, etc")
        synapse_channel = input()
    
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
    for max_file in max_files:
        print(os.path.basename(max_file)[:-4])
        color_files = glob.glob(os.path.basename(max_file)[:-4] + "*color.gif")
        
        if len(color_files) == 0:
            print("Couldn't find annotated file!")
    
        else:
            prints = []
            synapse_prints = []
            skeletons = []
        
            for color_file in color_files:
                some_prints, some_synapse_prints, some_skeletons = process_img(max_file, color_file)
                
                prints += some_prints
                synapse_prints += some_synapse_prints
                skeletons += some_skeletons
      
            #Save skeletons in "/processed/FOV#" folder
            for i in range(len(skeletons)):
                Image.fromarray(skeletons[i]).save(skel_dir +  "/" + os.path.basename(max_file)[:-4] + "_skel_" + str(i+1) + ".gif")
                Image.fromarray(prints[i]).save(skel_dir +  "/" + os.path.basename(max_file)[:-4] + "_print_" + str(i+1) + ".gif")
                Image.fromarray(synapse_prints[i]).save(skel_dir +  "/" + os.path.basename(max_file)[:-4] + "_synpase_print_" + str(i+1) + ".gif")

            # Getting mrna prints 
            for mRNA_channel in mRNA_channels:
                print("Channel: " + mRNA_channel)
                end_index = max_file.index(map2_channel)
                print(os.path.basename(max_file)[4:end_index] + mRNA_channel + ".tif")
                raw_file = glob.glob(os.path.basename(max_file)[4:end_index] + mRNA_channel + ".tif")[0]
                dapi_file = glob.glob(os.path.basename(max_file)[4:end_index] + dapi_channel + ".tif")[0]
                  
                #Making new subfolder to store outlines each channel
                channel_dir = "Outlines_" + mRNA_channel
                if not os.path.exists(channel_dir):
                  os.makedirs(channel_dir) 
                  
                getOutlinesTxts(raw_file, dapi_file, channel_dir, prints)
            
            #Making new subfolder to store outlines each channel
                synapse_dir = "Outlines_Synapses"
                if not os.path.exists(synapse_dir):
                  os.makedirs(synapse_dir) 

            # Getting synapse outlines
            synapse_file = glob.glob(os.path.basename(max_file)[4:end_index] + synapse_channel + ".tif")[0]
            getOutlinesTxts(synapse_file, dapi_file, synapse_dir, synapse_prints)
    
    return

def process_img(file, color_file):

  #reading in colored image and regular image
  color_img = plt.imread(color_file)[:,:,0:3]
  img = plt.imread(file)
  
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


  skeletons = [None]*num_colors
  outlines = [None]*num_colors
  synapse_outlines = [None]*num_colors
  
  #extracting individual dendrites from image
  j = 0
  for color in colors:
    if j < num_colors: #don't want to look for colors which are not
      mask = np.where(np.all(color_img == color, axis=-1), 1, 0)
      dendrite = np.multiply(img, mask)
      outlines[j], synapse_outlines[j], skeletons [j] = getDendriteAndSkeleton(dendrite, img)
      

      #trimming skeleton
      skeleton = skeletons[j]
      fil = FilFinder2D(skeleton, distance=1 * u.pix, mask=skeleton)
      #fil.preprocess_image(flatten_percent=85)
      fil.create_mask(border_masking=True, verbose=False,
                      use_existing_mask=True)
      fil.medskel(verbose=False)
      fil.analyze_skeletons(branch_thresh=100* u.pix, skel_thresh=40 * u.pix, prune_criteria='length')
      #parameters can be adjusted
      skel = fil.skeleton_longpath
      
      #get rid of extra circles next to skeleton
      skel_components, num_components = measure.label(skel, return_num=True)
      max_component_size = 0
      for label in range(1, num_components+1):
        skel_component = (skel_components == label)
        size = len(np.where(skel_component == 1))
        if size >= max_component_size:
          skel = skel_component
      
      skel[skel == 1] = 255
      skeletons[j] = skel

      j+= 1
  
  return outlines, synapse_outlines, skeletons

def getDendriteAndSkeleton(dendrite, img):
  
  #initial thresholding using whole image
  thresh = threshold_otsu(img)
  binary = dendrite > thresh

  #initial dilation of binary image of selected dendrite
  radius = 5
  selem = disk(radius)
  outline = binary_dilation(binary, selem=selem)

  #Repeatedly dilate ouline until binary image of dendrite is completely connected
  outline_copy = np.copy(outline)
  labeled, num_components = label(outline_copy, background=0, return_num=True)

  while num_components > 1:
    for i in range(10):
      outline_copy = binary_dilation(outline_copy, selem=disk(4))
      outline_copy = binary_erosion(outline_copy, selem=disk(3))
    labeled, num_components = label(outline_copy, background=0, return_num=True)

  #Skeletonize dilate outline
  skeleton = skeletonize(outline_copy)

  #This step only changes image if gaps in binary dendrite
  #Extract pixels of dilated skeleton not in outline
  adding = binary_dilation(skeleton, selem=disk(5))
  adding = adding * (adding != outline)
  outline = adding + outline #Adding to outline to fill gaps

  #Repeatedly dilate and erode to smooth outline
  for i in range(10):
    outline = binary_dilation(outline, selem=disk(3))
    outline = binary_erosion(outline, selem=disk(3))
  
  radius = 8
  selem = disk(radius)
  synapse_outline = binary_dilation(outline, selem=selem)
  
  return outline, synapse_outline, skeleton


def getOutlinesTxts(channel_img, dapi_img, channel_dir, prints):

  
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

  for i in range(len(prints)):

    #getting counterclockwise list of pixels on outlines
    contours = measure.find_contours(prints[i], 0)
    contours = np.array(contours)
    if np.shape(contours)[0] >= 2:
        print("This outline is sketchy!")
        """
        if len(contours[0]) > len(contours[1]):
            contours = np.array([contours[0]])
        else:
            contours = np.array([contours[1]])
        """
        contours = np.array([find_max_list(contours)])
    elif np.shape(contours)[0] == 0:
        print("Colours are wrong in this image!")
        break
    contours = contours.reshape(contours.shape[1:])

    outline_txt.write("CELL_START\tCell_" + str(i+1) + "\n")

    outline_txt.write("X_POS")
    for i in range(len(contours)):
     if i % 10 == 0:
        outline_txt.write("\t" + str(int(contours[i][1])))
    outline_txt.write("\t")

    outline_txt.write("\nY_POS")
    for i in range(len(contours)):
      if i % 10 == 0:
        outline_txt.write("\t" + str(int(contours[i][0])))
    outline_txt.write("\t")

    outline_txt.write("\nZ_Pos\t\n")
    outline_txt.write("CELL_END\n")

  outline_txt.close()

  return


def find_max_list(list):
    list_len = [len(i) for i in list]
    return list[np.argmax(np.array(list_len))]

main()