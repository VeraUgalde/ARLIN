# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 17:23:52 2021

@author: RyanHuang
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import os
import random
import seaborn as sns
from matplotlib.ticker import PercentFormatter
import pandas as pd


"""Get input data"""

def getmRNA(lines):
    mRNAs = []
    for line in lines:
        words = line.split("\t")
        y_coord = float(words[1])
        x_coord = float(words[0])
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


"""Calculate Euclidean distance between two points."""

def dist(pt, sig):
  if len(pt) != 2: print("AHHH")
  if len(sig) != 2: print("NOOOO")
  return  math.sqrt((pt[0] - sig[0])**2 + (pt[1] - sig[1])**2)

""" mRNA Colocalization"""
    
# Find the closest distance from a mRNA to a set of other mRNAs
def sigToOthers(sig, points, self_coloc=False):
    if self_coloc:
        k = np.min(np.array([dist(pt, sig) for pt in points if np.all(pt != sig)]))
    else: 
        k = np.min(np.array([dist(pt, sig) for pt in points]))
    return k

# Obtain a list of closest mRNA distances to other mRNAs in the same dendrite
def mRNAtomRNAs(mRNAs, othermRNAs):
    mRNADist = []

    #check if we are colocalizing mRNA with itself
    same = False
    if len(mRNAs)==len(othermRNAs): 
         same = np.all(mRNAs == othermRNAs)
         
    if same == True and len(mRNAs) <= 1:
        mRNADist = []
    else:
        for i in range(len(mRNAs)):
            mRNADist.append(sigToOthers(mRNAs[i], othermRNAs, same))
            
    return mRNADist

# Use positions of both mRNAs to analyze their colocalization
def mRNAcolocalization(mRNA1Positions, mRNA2Positions):
    mRNA12DistList = []
    mRNA21DistList = []
    for i in range(len(mRNA1Positions)):
        mRNA12Dist = mRNAtomRNAs(mRNA1Positions[i], mRNA2Positions[i])
        mRNA12DistList += mRNA12Dist
        mRNA21Dist = mRNAtomRNAs(mRNA2Positions[i], mRNA1Positions[i])
        mRNA21DistList += mRNA21Dist
    return mRNA12DistList, mRNA21DistList

#Get list of positions for control
def ctrlPosGeneration(num, dendrite):
    ctrlPixels  = [random.choice(dendrite) for i in range(num)]
    ctrlnanomt = []
    for pix in ctrlPixels:
        x_min = int(math.ceil(pix[0]*107.5))
        x_max = int(math.ceil(pix[0]*107.5 + 107.5))
        y_min = int(math.ceil(pix[1]*107.5))
        y_max = int(math.ceil(pix[1]*107.5 + 107.5))
        x_coord = random.choice(range(x_min, x_max))
        y_coord = random.choice(range(y_min, y_max))
        ctrlnanomt.append((x_coord, y_coord))
    return ctrlnanomt

""" Data organization and visualization """

# Plot a histogram for mRNA Colocalization (include all data)
def histogram(distances, mRNA, othermRNA, incre=20, x_lim=300):
    plt.clf()
    
    fig, axes = plt.subplots(5, 2,  figsize=(20, 35))
    
    y_lim = 0 # initialize hieght of all subplots
    
    for i in range(5):
        for j in range(2):
            # calculate number of bins according to the desired increment
            groups = int(max(distances[i][j])//incre)
    
            n, bins, patches = axes[i][j].hist(distances[i][j], groups, weights= np.ones(len(distances[i][j])) / len(distances[i][j]), edgecolor='blue')
            height = max(n)
            if height > y_lim: y_lim = height 
            axes[i][j].set_xlabel('Distance to Closest ' + othermRNA[i][j] + ' mRNA')
            axes[i][j].set_ylabel(mRNA[i][j] + ' mRNA Percentage')
            axes[i][j].set_title(mRNA[i][j] + ' mRNA Distribution to Closest ' + othermRNA[i][j])


    #align y axes of all subplots
    for i in range(5):
        for j in range(2):
            axes[i][j].yaxis.set_major_formatter(PercentFormatter(1))
            axes[i][j].set_ylim(0, y_lim + 0.005)

    plt.savefig("mRNAColocHistogram.png")

    
# Organize and output the mRNA Colocalization data to a table
def outputTable(distances, incre=20):
    
    totalmRNAcount = len(distances)
    groups = int(math.ceil(np.max(distances)/incre))
    mod_distances = [int(dist//20) for dist in distances]
    mod_distances = np.array(mod_distances)
    
    cols = ["Distance (nm)", "mRNA Count", "Percentage"]
    dataframe = pd.DataFrame(columns=cols)
    
    x_min = 0
    x_max = incre
    
    for j in range(groups):
        distance_group = str(x_min) + '-' + str(x_max)
        count = np.count_nonzero(mod_distances == j)
        dataframe2 = pd.DataFrame([[distance_group, count, count/totalmRNAcount * 100]], columns=cols)
        dataframe = dataframe.append(dataframe2, ignore_index=True)
        x_min += incre
        x_max += incre
    
    return dataframe

# Plot a barplot for mRNA Colocalization (select for desired range)
def makePlot(data, titles, maxD = 300, incre = 20):
    
    plt.clf()
    
    fig, axes = plt.subplots(5, 2,  figsize=(30, 45))
     
    y_lim = 0 # initialize hieght of all subplots
    for i in range(5):
        for j in range(2):
            df = data[i][j].head(int(maxD//incre))
            height = df["Percentage"].max()
            if height > y_lim: y_lim = height
            sns.barplot(x="Distance (nm)", y="Percentage", data=df, ax=axes[i][j]).set_title(titles[i][j])
            print(titles[i][j])

    #align y axes of all subplots
    for i in range(5):
        for j in range(2):
            axes[i][j].set_ylim(0, y_lim + 0.5)
            
    plt.savefig("ColocalizationBarplots.png")
    

""" Main Program """
 
def main():
    
    # Getting the image folder
    print("\n Full path to experiment folder: ")
    directory = input()
    
    while not os.path.isdir(directory):
      print("\n Incorrect path, please try again:")
      print("\n Full path to folder: ")
      directory = input()
    
    # Getting mRNA channels and plotting parameters
    print("What are the channels you want to analyze for colocalization?")
    print("Enter the first channel now:")
    channel1 = input()
    print("Enter the second channel now:")
    channel2 = input()
    mRNA1PosLists = []
    mRNA2PosLists = []
    ctrl1PosLists = [[] for k in range(50)]
    ctrl2PosLists = [[] for k in range(50)]
    
    print("Enter the mRNA corresponding to " + channel1 + " for colocalization?")
    mRNA1 = input()
    print("Enter the mRNA corresponding to " + channel2 + " for colocalization?")
    mRNA2 = input()
    
    print("Enter the maximum distance you want to plot: ")
    max_dist = int(input())
    print("Enter the desired increment for distance on the plots: ")
    dist_increment = int(input())
    

    # Getting MAP2 channel and check for skeletons
    map2_channel = " " 
    print("\n Name of MAP2 channel: ")
    print("\n ex: CY7, c7, DAPI, etc")
    map2_channel = input()
    while not glob.glob("*/SkeletonsAndPrints/MAX_*" + map2_channel + "*skel_*" + ".gif"):
        print("\n Incorrect channel name. Try again: ")
        print("\n Name of MAP2 channel: ")
        print("\n ex: CY7, c7, DAPI, etc")
        map2_channel = input()
        
    
    
    # Process dendrites by looping through both channels and each field of view
    files = glob.glob(directory + "/Results_" + channel1 + "/*xy*outline_spots*" + ".txt")

    div_treats = []
    for file in files:
        endindex = os.path.basename(file).index("__outline_spots")
        file_name = os.path.basename(file)[:endindex]
        print(file_name)
        div = file_name.split('_')[1][3:]
        treatment = file_name.split('_')[2]
        if (div,treatment) not in div_treats:
            div_treats.append((div,treatment))

    for div, treatment in div_treats:
    
        for file in files:
            endindex = os.path.basename(file).index("__outline_spots")
            file_name = os.path.basename(file)[:endindex]
            print(file_name)
            f_div = file_name.split('_')[1][3:]
            f_treatment = file_name.split('_')[2]
            if f_div == div and f_treatment == treatment:
                channel_index = file_name.index(channel1)
                matched_file = glob.glob(directory + "/Results_" + channel2 + "/" + file_name[:channel_index] + channel2 + "*outline_spots*" + ".txt")[0]
                
                dendrites1 = inputmRNA(file)
                dendrites2 = inputmRNA(matched_file)
                
                j = 1
                
                print_files = glob.glob(directory + "/SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "*print_*" + ".gif")
                print(len(print_files))
                print(len(dendrites1))
                
                for i in range(len(dendrites1)):
                    
                    if len(dendrites1[i]) != 0 and len(dendrites2[i]) != 0:
                    
                        #Getting print file and generating skeleton points
                        print_file = glob.glob(directory + "/SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "*print_" + str(j) + ".gif")[0]
                        dendrite_print = plt.imread(print_file)
                        print_points = np.argwhere(dendrite_print == 255)
                        
                        while len(print_points) == 0:
                            j = j+1
                            print_file = glob.glob(directory + "/SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "*print_" + str(j) + ".gif")[0]
                            dendrite_print = plt.imread(print_file)
                            print_points = np.argwhere(dendrite_print == 255)
                    
                        # For each dendrite, add mRNA1/mRNA2 data to the lists
                        mRNA1PosLists.append(np.array(dendrites1[i]))
                        for k in range(50):
                            ctrl1PosLists[k].append(ctrlPosGeneration(len(dendrites1[i]), print_points))
                            ctrl2PosLists[k].append(ctrlPosGeneration(len(dendrites2[i]), print_points))
                        mRNA2PosLists.append(np.array(dendrites2[i]))
                        
                    
                    j = j+1
        
        #print([len(mRNA1PosLists[i]) for i in range(len(mRNA1PosLists))])
        #print([len(mRNA2PosLists[i]) for i in range(len(mRNA2PosLists))])
        
        # Process the data from all dendrites
        distances = []
        temp1 = []
        temp2 = []
        distances.append(mRNAcolocalization(mRNA1PosLists, mRNA2PosLists))
        for i in range(50):
            dist12, dist21 = mRNAcolocalization(ctrl1PosLists[i], mRNA2PosLists)
            temp1 += dist12
            temp2 += dist21
        distances.append([temp1, temp2])
        
        temp1 = []
        temp2 = []
        for i in range(50):
            dist12, dist21 = mRNAcolocalization(mRNA1PosLists, ctrl2PosLists[i])
            temp1 += dist12
            temp2 += dist21
        distances.append([temp1, temp2])
        
        temp1 = []
        temp2 = []
        for i in range(50):
            for j in range(50):
                dist12, dist21 = mRNAcolocalization(ctrl1PosLists[i], ctrl2PosLists[j])
                temp1 += dist12
                temp2 += dist21
        distances.append([temp1, temp2])
        
        mRNA11Dists, unimportant = mRNAcolocalization(mRNA1PosLists, mRNA1PosLists)
        mRNA22Dists, unimportant = mRNAcolocalization(mRNA2PosLists, mRNA2PosLists)
        distances.append([mRNA11Dists, mRNA22Dists])
        
        mRNAs = [[mRNA1, mRNA2], ["Ctl-" + mRNA1, mRNA2], [mRNA1, "Ctl-" + mRNA2], ["Ctl-" + mRNA1, "Ctl-" + mRNA2], [mRNA1, mRNA2]]
        othermRNAs = [[mRNA2, mRNA1], [mRNA2, "Ctl-" + mRNA1], ["Ctl-" + mRNA2, mRNA1], ["Ctl-" + mRNA2, "Ctl-" + mRNA1], [mRNA1, mRNA2]]
        
        
        # Generate histograms
        histogram(distances, mRNAs, othermRNAs, dist_increment)

        tables = [[None]*2 for k in range(5)]
        titles = [[None]*2 for k in range(5)]
        with pd.ExcelWriter("Colocalization-" + div + "-" + treatment + ".xlsx") as writer:
            for i in range(5):
                for j in range(2):
                    table = outputTable(distances[i][j], dist_increment)
                    tables[i][j] = table
                    title = mRNAs[i][j] + " to closest " + othermRNAs[i][j]
                    titles[i][j] = title
                    table.to_excel(writer, sheet_name=title)
        
        makePlot(tables, titles, max_dist, dist_increment) # make bar plots


    
main()