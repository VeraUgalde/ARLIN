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
from collections import defaultdict

"""Calculate Euclidean distance between two points."""


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

        
    density_dict = defaultdict(list())
    
    os.chdir(directory)
    

    for channel in mRNA_channels:
        files = glob.glob("Results_" + channel + "/*XY*outline_spots*" + ".txt")
        
        
        for file in files:
            endindex = os.path.basename(file).index("__outline_spots")
            file_name = os.path.basename(file)[:endindex]
            print(file_name)
            div = file_name.split('_')[1][3:]
            treatment = file_name.split('_')[2]
            sheet_name = "DIV" + div + " " + treatment
            somas = inputmRNA(file)
            print(len(somas))
            
            for i in range(len(somas)):
                
                #Getting skeleton file and generating skeleton points
                channel_index = file_name.index(channel)
                print_files = glob.glob("SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "*print_" + str(i+1) + ".gif")
                nuc_print_files = glob.glob("SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "*nuc_print_" + str(i+1) + ".gif")
                if len(print_files) == 0:
                    print("Missing print for " + file_name)
                    break

                
                print_file = print_files[0]
                soma_print = plt.imread(print_file)
                print_points = np.argwhere(soma_print == 255)
                signals = np.array(somas[i])

                if len(nuc_print_files) != 0:
                    nuc_print_file = nuc_print_files[0]
                    nuc_print = plt.imread(nuc_print_file)
                    nuc_points = np.argwhere(nuc_print == 255)
                    if len(nuc_points) == 0:
                        nuc_area = ""
                    else:
                        nuc_area = len(nuc_points)
                        

                if len(print_points) != 0:
                
                    new_density = [channel, len(signals), len(print_points), len(signals)/len(print_points), nuc_area]
                    density_dict[sheet_name].append(new_density)
    
        i = 0
    
    with pd.ExcelWriter("SomaDensities.xlsx") as writer: #change name later
        for sheetname, densities in density_dict.items():
            densities = pd.DataFrame(densities, columns=["Channel", "Num mrna", "Area", "Density", "Nucleus Area"])
            densities.to_excel(writer, sheet_name=sheetname, columns=["Channel", "Num mrna", "Area", "Density", "Nucleus Area"], index=False)


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


main()