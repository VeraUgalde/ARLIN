#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import math
import matplotlib.pyplot as plt
import glob
import os
import random
import seaborn as sns
from matplotlib.ticker import PercentFormatter
import pandas as pd
from collections import Counter
from collections import defaultdict


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
  return  math.sqrt((pt[0] - sig[0])**2 + (pt[1] - sig[1])**2)

""" mRNA Colocalization"""

def distanceCount(synap_pts, mrna_pts, dist_thresh):
    counts = []
    for synap in synap_pts:
        count = np.count_nonzero([dist(synap, mrna) < dist_thresh for mrna in mrna_pts])
        counts.append(count)
    return counts


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
            axes[i][j].set_ylim(0, y_lim)
            

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
def makeGraph(mrna_counts, ctrl_counts, mrna_name):
    
    print(len(mrna_counts), len(ctrl_counts))
    xlim = max(max(mrna_counts), max(ctrl_counts))
    print(xlim)

    plt.clf()
    
    fig, axes = plt.subplots(1, 2,  figsize=(30, 15))
    axes[0].set_xlim(0, xlim)
    axes[1].set_xlim(0, xlim)
    

    n_mrna, bins, patches = axes[0].hist(mrna_counts, bins=5, weights=np.ones(len(mrna_counts)) / len(mrna_counts))
    axes[0].set_title(mrna_name + " mRNA Colocalization with Synapse")
    axes[0].set_xlabel("Number of colocalized mRNA")
    axes[0].set_ylabel("Percent of Synapses")

    n_ctrl, bins, patches = axes[1].hist(ctrl_counts, bins=5, weights=np.ones(len(ctrl_counts)) / len(ctrl_counts))
    axes[1].set_title("Random Simulated Colocalization with Synapse")
    axes[1].set_xlabel("Number of colocalized Simulated Points")
    axes[1].set_ylabel("Percent of Synapses")   

    print(len(n_mrna), len(n_ctrl))

    #Setting x/y axis range on subplots
    ylim = max([max(n_mrna), max(n_ctrl)])
    axes[0].set_ylim(0, ylim + 0.005)
    axes[1].set_ylim(0, ylim +0.005)

    #Setting yaxis to be percent
    axes[0].yaxis.set_major_formatter(PercentFormatter(1))
    axes[1].yaxis.set_major_formatter(PercentFormatter(1))

    plt.savefig(mrna_name + "_Synapse_Colocalization.png")
    

""" Main Program """
 
def main():
    
    # Getting the image folder
    print("\n Full path to experiment folder: ")
    directory = input()
    
    while not os.path.isdir(directory):
      print("\n Incorrect path, please try again:")
      print("\n Full path to folder: ")
      directory = input()
    
    # Getting synapse/mrna channels and plotting parameters
    print("What is the synapse channel? \n")
    synap_channel = input()

    print("How many mRNA channels are you looking at for colocalization? \n")
    num = int(input())
    mrnas = []
    for i in range(num):
        print("Name of the mRNA channel: \n")
        channel = input()
        print("Name of mRNA corresponding to " + channel + "\n")
        name = input()
        mrnas.append((channel, name))
    

    print("Enter the distance threshold for colocalization ")
    dist_thresh = int(input())
    

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
        
    
    synap_files = glob.glob(directory + "/Results_" + synap_channel + "/*XY*outline_spots*" + ".txt")

    
    sheets = {}

    synapse_dendrites_list = []
    prints_list = []
    for synap_file in synap_files:

        #Getting skeleton file and generating skeleton points
        endindex = os.path.basename(synap_file).index("__outline_spots")
        file_name = os.path.basename(synap_file)[:endindex]
        print(file_name)
        div = file_name.split('_')[1][3:]
        treatment = file_name.split('_')[2]
        sheet_name = "DIV" + div + " " + treatment
        channel_index = file_name.index(synap_channel)

        #Getting synapse coordinates
        synapse_dendrites = inputmRNA(synap_file)

        synapse_dendrites_list.append(synapse_dendrites)
        
        prints = []
        print_files = glob.glob(directory + "/SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "_print_*.gif")
        print(len(print_files))
        for i in range(len(print_files)):
            #Get dendrite print
            print_file = glob.glob(directory + "/SkeletonsAndPrints" + "/MAX_" + file_name[:channel_index] + map2_channel + "_print_" + str(i+1) + ".gif")[0]
            dendrite_print = plt.imread(print_file)
            print_points = np.argwhere(dendrite_print == 255)
            if len(print_points) != 0:
                print(i)
                prints.append(print_points)

        prints_list.append(prints)

    dfs = []
    for (mrna_channel, mrna_name) in mrnas:
        mrna_counts = defaultdict(Counter)
        ctrl_counts = defaultdict(Counter)
        
        print("second loop")
        for i in range(len(synap_files)):
            endindex = os.path.basename(synap_files[i]).index("__outline_spots")
            file_name = os.path.basename(synap_files[i])[:endindex]
            div = file_name.split('_')[1][3:]
            treatment = file_name.split('_')[2]
            sheet_name = "DIV" + div + " " + treatment
            print(file_name)
            matched_file = glob.glob(directory + "/Results_" + mrna_channel + "/" + file_name.replace(synap_channel, mrna_channel) + "*outline_spots*" + ".txt")[0]
            mrna_dendrites = inputmRNA(matched_file)
            print(len(mrna_dendrites))
            print(len(prints_list[i]))

            for j in range(len(mrna_dendrites)):
                #Generate control signal
                
                ctrl_runs = Counter([])
                for k in range(50):
                    ctrl_points = ctrlPosGeneration(len(mrna_dendrites[j]), prints_list[i][j])
                    ctrl_counts_run = distanceCount(synapse_dendrites_list[i][j], ctrl_points, dist_thresh)
                    ctrl_runs.update(ctrl_counts_run)
                    
                for k in ctrl_runs: ctrl_runs[k] /= 50
                    
                    
                #Find distances between points and count those less than threshold
                ctrl_counts[sheet_name].update(ctrl_runs)
                mrna_counts[sheet_name].update(distanceCount(synapse_dendrites_list[i][j], mrna_dendrites[j], dist_thresh))
        

        max_num = max(list(mrna_counts.keys()) + list(ctrl_counts.keys()))
        mrna_data = [0 if i not in mrna_counts.keys() else mrna_counts[i] for i in range(max_num)]
        ctrl_data = [0 if i not in ctrl_counts.keys() else ctrl_counts[i] for i in range(max_num)]
        counts = [i for i in range(max_num)]
        data = np.array([counts, mrna_data, ctrl_data]).T
        df = pd.DataFrame(data, columns=["Num of mRNA", "Real " + mrna_name, "Sim " + mrna_name])
        dfs.append(df)
       # makeGraph(mrna_counts, ctrl_counts, mrna_name)
    df = pd.concat(dfs, axis=1)
    sheets[sheet_name] = df

    with pd.ExcelWriter("synapse_coloc.xlsx") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)

main()