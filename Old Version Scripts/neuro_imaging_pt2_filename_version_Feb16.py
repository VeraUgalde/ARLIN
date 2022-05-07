# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 14:23:27 2022

@author: Administrator
"""


import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import os
import xlsxwriter
import pandas as pd
import random

"""Calculate Euclidean distance between two points."""

def dist(pt, sig):
  return  math.sqrt((pt[0] - sig[0])**2 + (pt[1] - sig[1])**2)

def main():
    print("\n Full path to experiment folder: ")
    print("\n 21-01-29 neurons Ctrl-HS-MG132_eEF1A2-HSPA1A-MAP2-PSD95")
    directory = input()
    
    while not os.path.isdir(directory):
      print("\n Incorrect path, please try again:")
      print("\n Full path to folder: ")
      directory = input()
    
    #Getting mRNA channels
    mRNA_channel = " "
    mRNA_channels = []
    while mRNA_channel != "0":
      print("\n Name of mRNA channel (if done press 0): ")
      print("\n ex: CY5, c5, CY3, c3 etc")
      mRNA_channel = input()
      while not glob.glob("*/*XY*" + mRNA_channel + ".tif") and mRNA_channel != "0": 
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
    while not glob.glob("*/*XY*" + map2_channel + ".tif"):
        print("\n Incorrect channel name. Try again: ")
        print("\n Name of MAP2 channel: ")
        print("\n ex: CY7, c7, DAPI, etc")
        map2_channel = input()

        
    density_dict = {}
    
    os.chdir(directory)
    
    cols = ["mRNA Channel", "DIV", "Treatment", "Soma Distance", "Dendrite", "mRNA Count"]
    range_dict = {0 : "0-25", 1 : "25-50", 2 : "50-75", 3 : "75-100", 4 : "100-125", 5 : "125-150"}
    
    with pd.ExcelWriter("Data_sb_new.xlsx") as writer: #change name later   
        dataframes = []
        sheet_names = []
    
        for channel in mRNA_channels:
            files = glob.glob("Results_" + channel + "/*XY*outline_spots*" + ".txt")
            
            dendrite_counts = []
            for p in sheet_names:
                dendrite_counts.append(0)
            
            for file in files:
                endindex = os.path.basename(file).index("__outline_spots")
                file_name = os.path.basename(file)[:endindex]
                print(file_name)
                div = file_name.split('_')[1][3:]
                treatment = file_name.split('_')[2]
                sheet_name = "DIV" + div + " " + treatment
                if sheet_name not in sheet_names:
                    dataframe = pd.DataFrame(columns=cols)
                    dataframes.append(dataframe)
                    dataframe_index = len(dataframes) - 1
                    sheet_names.append(sheet_name)
                    dendrite_count = 0
                    dendrite_counts.append(dendrite_count)
                    density_dict[sheet_name] = []
                else:
                    dataframe_index = sheet_names.index(sheet_name)
                    dataframe = dataframes[dataframe_index]
                    dendrite_count = dendrite_counts[dataframe_index]
                
                dendrites = inputmRNA(file)
                print(len(dendrites))
                
                for i in range(len(dendrites)):
                    
                    
                    #Getting skeleton file and generating skeleton points
                    channel_index = file_name.index(channel)
                    skeleton_files = glob.glob("SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "*skel_" + str(i+1) + ".gif")
                    print_files = glob.glob("SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "*print_" + str(i+1) + ".gif")
                    if len(skeleton_files) == 0 or len(print_files) == 0:
                        print("Missing print/skel for " + file_name)
                        break
                    
                    print_file = print_files[0]
                    skeleton_file = skeleton_files[0]
                    if len(np.array(plt.imread(skeleton_file)).shape) == 3:
                        skeleton = np.squeeze(plt.imread(skeleton_file)[:,:,0])
                    else:
                        skeleton = plt.imread(skeleton_file)
                    dendrite_print = plt.imread(print_file)
                    skel_points = np.argwhere(skeleton == 255)
                    print_points = np.argwhere(dendrite_print == 255)
                    
                  
                    
                    signals = np.array(dendrites[i])
                    if len(print_points) != 0:
                    
                        new_density = [channel, len(signals), len(print_points), len(signals)/len(print_points)]
                        density_dict[sheet_name].append(new_density)
              
                        endpoints = findEndPoints(skeleton, skel_points)
    
    
                        if len(endpoints) != 2:
                            print("Skeleton from "  + skeleton_file + " doesn't have 2 endpoints")
                        elif len(print_points) == 0:
                            print("Print from " + print_file + " is empty")
                        else:
                            #dist_array is list of 1s and root-2s
                            ordered_skel, dist_array = orderSkelPoints(skeleton, skel_points, findSoma(signals, endpoints))
                            mapped_signals = mapping(ordered_skel, signals)
                        
                            sig_distances = []
              
                            for sig in mapped_signals:
                                k = np.argwhere(np.all(ordered_skel == sig, axis=1))[0][0]
                                sig_distances.append((np.sum(dist_array[0:k+1]) * 0.1075) // 25)
                          
                            dendrite_length = np.sum(dist_array) * 0.1075
                            print(dendrite_length)
              
                            sig_distances = np.array(sig_distances)
                            for j in range(6):
                                if dendrite_length > (j)*25:
                                    count = np.count_nonzero(sig_distances == j)
                                    dataframe2 = pd.DataFrame([[channel, int(div), treatment, range_dict[j], dendrite_count + 1, count]], columns=cols)
                                    dataframe = dataframe.append(dataframe2, ignore_index=True)
                                    dataframes[dataframe_index] = dataframe
                        
                            dendrite_count = dendrite_count + 1
                            dendrite_counts[dataframe_index] = dendrite_count
        
        i = 0  
        for dataframe in dataframes:            
            dataframe.to_excel(writer, sheet_name=sheet_names[i], columns=cols, index=False)
            i = i + 1
    
    with pd.ExcelWriter("Densities_new.xlsx") as writer: #change name later
        for sheetname, densities in density_dict.items():
            densities = pd.DataFrame(densities, columns=["Channel", "Num mrna", "Area", "Density"])
            densities.to_excel(writer, sheet_name=sheetname, columns=["Channel", "Num mrna", "Area", "Density"], index=False)
       


def findEndPoints(skel, points):
  endpoints = []
  for point in points:
    if np.sum((skel[(point[0] - 1):(point[0] + 2), (point[1] - 1):(point[1] + 2)])) == 510:
      endpoints.append(point)
  return endpoints

"""Find the endpoint that corresponds to the soma."""

def findSoma(mRNAs, endpoints):
    if len(endpoints) == 2:
        point1 = endpoints[0]
        point2 = endpoints[1]
        totaldistance1 = 0
        totaldistance2 = 0
        for point in mRNAs:
            totaldistance1 = totaldistance1 + dist (point, point1)
            totaldistance2 = totaldistance2 + dist (point, point2)
        somaend = point1
        if totaldistance2 < totaldistance1:
            somaend = point2
    elif len(endpoints) == 1:
        somaend = endpoints[0]
    else:
        print("This dendrite is problematic!")
        somaend = None
        
    return somaend

"""Create an ordered list of skeleton pixels."""

def orderSkelPoints(skel, skel_points, endpoint):
  dist_array = np.zeros(len(skel_points))
  skel_array = np.zeros(skel_points.shape)
  prev_coord = endpoint
  cur_coord = endpoint

  for i in range(len(skel_points)):
    skel_array[i] = cur_coord
    dist_array[i] = dist(prev_coord, cur_coord)

    pixels = np.argwhere(skel[cur_coord[0] - 1:cur_coord[0] + 2, cur_coord[1] - 1:cur_coord[1] + 2] == 255)
    pixels = [(cur_coord - (1,1) + pix) for pix in pixels]

    cur_coord_copy = cur_coord
    prev_coord_copy = prev_coord
    
    for pix in pixels:
      #print("pix " + str(pix))
      if np.any(pix != cur_coord_copy) and np.any(pix != prev_coord_copy):
        # and skel[pix[0], pix[1]] == 255
        prev_coord = cur_coord
        cur_coord = pix
        
  
  return skel_array, dist_array

"""Map single mRNA signal to its corresponding point on skeleton."""

def mapSig(sig, points):
  k = np.argmin(np.array([dist(pt, sig) for pt in points]))
  return points[k]

"""Find the distance of mRNA to a closest synapse."""

def sigToSynap(sig, points):
  k = np.min(np.array([dist(pt, sig) for pt in points]))
  return k

"""Create a list with all mRNA signals mapped onto their corresponding points on skeleton."""

def mapping(skel_points, signals):
  mapped_signals = np.zeros(signals.shape)
  for i in range(len(signals)):
    mapped_signals[i,:] = mapSig(signals[i], skel_points)
  return mapped_signals

"""Output raw data to a textfile."""

def outputRawData(dendrite_length, sig_distances):
  
    raw_txt = open("raw_data.txt", "a")

    raw_txt.write(str(dendrite_length) + "\n")
    
    for sig_distance in sig_distances:
        raw_txt.write(str(sig_distance) + "\n")

    raw_txt.close()
    return


"""Process a field of view: mRNA distance to soma"""

def processFOV_soma(num_FOV, dendrites):
    
    dendr_num = 1
    
    for dendrite in dendrites:
        
        skel_file = glob.glob("MAX*XY" + str(num_FOV) + "*skel_" + str(dendr_num) + ".gif")
        skel = plt.imread(skel_file[0])
        skel_points = np.argwhere(skel == 255)
        

        signals = np.array(dendrite)
        print(signals)

        endpoints = findEndPoints(skel, skel_points)
        if len(endpoints) != 2:
            print("\nDendrite ")

        mapped_signals = map(skel_points, signals)
        print(mapped_signals)

        ordered_skel, dist_array = orderSkelPoints(skel, skel_points, findSoma(signals, endpoints))
        print(ordered_skel)


        sig_distances = []


        for sig in mapped_signals:
          k = np.argwhere(np.all(ordered_skel == sig, axis=1))[0][0]
          sig_distances.append(np.sum(dist_array[0:k+1]))

        dendrite_length = np.sum(dist_array)

        print(dendrite_length)
        print(sig_distances)

        return (dendrite_length, sig_distances)
        
        dendr_num = dendr_num + 1

"""Main program"""

def getmRNA(lines):
    mRNAs = []
    for line in lines:
        words = line.split("\t")
        y_coord = float(words[1])/107.5
        x_coord = float(words[0])/107.5
        mRNAs.append((x_coord, y_coord))
    return mRNAs

def inputmRNA(filename):
    file = open(filename, "r")

    startLines = []
    endLines = []
    dendrite_lines = []
    currentLine = 1

    for line in file:
        if line == "SPOTS_START" + "\n":
            startLines.append(currentLine)
        elif line == "SPOTS_END" + "\n":
            endLines.append(currentLine)
        elif line == "CELL_END" + "\n":
            dendrite_lines.append(currentLine)
        currentLine += 1

    file.close()
    
    
    file = open(filename, "r")
    lineLists = list(file)
    file.close()
    
    dendrites = []
    dendrite_with_dot = 0
    for dendrite in dendrite_lines:
        if dendrite + 1 in startLines:
            dendrites.append(getmRNA(lineLists[startLines[dendrite_with_dot] + 1 : endLines[dendrite_with_dot] - 1]))
            dendrite_with_dot += 1
        else:
            dendrites.append([])
        
    
    return dendrites



# mRNA Colocalization
    
# Find the closest distance from a mRNA to a set of other mRNAs
def sigToOthers(sig, points):
    k = np.min(np.array([dist(pt, sig) for pt in points]))
    return k

# Obtain a list of closest mRNA distances to other mRNAs in the same dendrite
def mRNAtomRNAs(mRNAs, othermRNAs):
    mRNADist = []
    for i in range(len(mRNAs)):
        mRNADist.append(sigToSynap(mRNAs[i], othermRNAs))
    return mRNADist



main()