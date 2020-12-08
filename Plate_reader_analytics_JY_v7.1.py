###############################################################################
# Plate_reader_analytics.py v7.1
# 
# Johannes Yayo
# December 8, 2020
# Script to calculate growth rate from optical density (OD) data generated
# in a 96-well plate format
#
# Instruction on use are found in the README-file
###############################################################################

#! /usr/bin/env python

# perform required imports
import sys
import os
import numpy as np
import pandas as pd
from bisect import bisect_right
import matplotlib.pyplot as plt
from scipy.stats import linregress
from timeit import default_timer as timer


class SlopePoint:
    '''
    This class will hold information on a slope
    slope   - value of the slope
    start   - start index for the slope
    end     - end index for the slope
    r_value - the statistic parameter R2
    '''
    def __init__(self, slope, start, end, r_value):
        self.slope = slope
        self.start = start
        self.end = end
        self.r_value = r_value

def get_index(serie, value):
    '''
    Returns the index of a specified value from a series
    '''
    for a in serie:
        if a >= value:
            return (serie[serie==a].index[0])
    return 'NaN'

def getIndicesFromTable(searchTerms, table):
    '''
    Finds the corresponding 96-well plate indices (row+column)
    searchTerms - the terms (string) of interest 
    table       - a 2D list with 8 x 12 elements of strings
    
    Returns a list with indices in the format row+column, with 
    row in letters (A-H) and column in numbers (1-12)
    '''
    sumOfIndices = []
    rows = ['A','B','C','D','E','F','G','H']
    
    for searchTerm in searchTerms:
        indexList = []
        for i in range(8):
            for j in range(12):
                if table[i][j] == searchTerm:
                    indexList.append(rows[i]+str(j+1))
        sumOfIndices.append(indexList)
        
    return sumOfIndices

def intersection(list1, list2):
    '''
    Finds the intersection between two lists, i.e. elements shared. 
    list1,2     - two lists of strings
    '''
    return [value for value in list1 if value in list2]

def superLabel(f, axis,label,label_prop=None, labelpad=7, ha='center',va='center'):
    ''' Add super ylabel or xlabel to the figure
    Similar to matplotlib.suptitle
    axis       - string: "x" or "y"
    label      - string
    label_prop - keyword dictionary for Text
    labelpad   - padding from the axis (default: 5)
    ha         - horizontal alignment (default: "center")
    va         - vertical alignment (default: "center")
    '''
    xmin = []
    ymin = []
    for ax in f.axes:
        xmin.append(ax.get_position().xmin)
        ymin.append(ax.get_position().ymin)
    xmin,ymin = min(xmin),min(ymin)
    dpi = f.dpi
    if axis.lower() == "y":
        rotation=90.
        x = xmin-float(labelpad-2)/dpi
        y = 0.5
    elif axis.lower() == 'x':
        rotation = 0.
        x = 0.5
        y = ymin - float(labelpad)/dpi
    else:
        raise Exception("Unexpected axis: x or y")
    if label_prop is None:
        label_prop = dict()
    f.text(x,y,label,rotation=rotation,
               transform=f.transFigure,
               ha=ha,va=va, fontsize=12,
               **label_prop)
    
def plotStrain(name, data, stats, pos, media, fpath, plotLN = False, draw_slopePoints = False, draw_lag_time = False, print_stats = False):
    '''
    Plots one grid with each column representing the media tested and each row
    represents different wells
    
    name                - strain name 
    data                - DataFrame with plot data incl time_h
    pos                 - 2D list with indexes of the strain on each media type
    media               - 1D list with name(s) of media types to be plotted
    stats               - DataFrame with growth statistics per well
    fpath               - file path to save fig
    plotLN              - true if axis should be adjusted for LN values
    draw_slopePoints    - plots the exponential phase in a different color
    draw_lag_time       - draw vertical line at t = lag_time
    print_stats         - prints mu-, R2- and lag time-values in the plot
    '''
    
    # Find number of columns and rows for the grid
    cols = len(media)
    rows = 0
    for i in range(cols):
        if len(pos[i]) > rows:
            rows = len(pos[i])

    # Create figure and set plotting details
    fig, ax = plt.subplots(nrows=rows,ncols=cols, sharex=True, figsize=(8,1.6*rows))
    
    markerSize = 1
    xMin = 0
    xMax = 80
    yMin = 0
    yMax = 1.5
    xTicks = np.linspace(0, xMax, 9)
    yTicks = np.linspace(yMin, yMax, 4)
    if plotLN:
        yMin = -5.5
        yMax = 0.5
        yTicks = np.linspace(yMin, yMax, 5)
    
    # Make each plot in the grid
    for i in range(cols):
        
        if (len(pos[i]) != 0): # no plotting if no wells are listed
            for j in range(rows):

                #extract column name (letter+number), i.e. well on 96-well plate
                well = pos[i][j]
            
                # Printing media name above first plot in column
                title = ''
                if j == 0:  
                    title = media[i]
                
                if len(media)==1:
                    plot(ax[j], data["time_h"], data[well], title, stats[well], markerSize, xMin,xMax,yMin,yMax,xTicks,yTicks,draw_slopePoints,draw_lag_time,print_stats)
                else:
                    plot(ax[j, i], data["time_h"], data[well], title, stats[well], markerSize, xMin,xMax,yMin,yMax,xTicks,yTicks,draw_slopePoints,draw_lag_time,print_stats)
    
    # Add a super x- and y-label for the axes
    lnText = "" 
    yText = "OD600"
    if plotLN:
        lnText = "- LN"
        yText = "LN(OD600)"
    
    fig.suptitle("Growth curves for "+name+lnText, fontsize=20)
    superLabel(fig, "y", yText)
    superLabel(fig, "x", "time (h)")
    fig.savefig(fpath)
    
def plot(ax,x,y,title,stats,markerSize,xMin,xMax,yMin,yMax,xTicks,yTicks,draw_slopePoints,draw_lag_time,print_stats):
    '''
    Plotting function with possibility to draw information on the plot
    ax                  - axis to plot on
    x,y                 - list with x and y values
    title               - string to print on top of plot, if not empty
    stats               - list with growth statistics (well, growth rate, R2, 
                          max ODlag, start pos, end pos, range)
    markerSize          - size of the marker for scatter plot
    xMin, xMax, xTicks  - x-axis settings
    yMin, yMax, yTicks  - y-axis settings
    draw_slopePoints    - plots the exponential phase in a different color
    draw_lag_time       - draw vertical line at t = lag_time
    print_stats         - prints mu-, R2- and lag time-values in the plot
    '''
    
    ax.scatter(x, y, s=markerSize, c='blue')                
    ax.set_xlim((xMin,xMax))
    ax.set_ylim((yMin, yMax))
    ax.set_xticks(xTicks)
    ax.set_yticks(yTicks)
    ax.text(0.65, 0.87, stats[0], fontsize=8, transform=ax.transAxes, color='red')  
    # Drawing the exponential phase
    if draw_slopePoints:  
        start = stats.iloc[5]
        end = stats.iloc[6]
        ax.scatter(x.iloc[start:end], y.iloc[start:end],
                                s=markerSize, c='red')
    # Draw a line representing lag_time
    if draw_lag_time:  
        lag_time = stats.iloc[4] 
        if (lag_time >= 0):
            ax.axvline(x=lag_time, dashes = (3,1))
    # Prints mu-, R2- and lag time-values in the plot
    if print_stats: 
        ax.text(0.65, 0.74, "$\mu$="+str(round(stats.iloc[1], 3)), fontsize=8, transform=ax.transAxes, color='red')  # print mu
        ax.text(0.65, 0.61, "R2=" +str(round(stats.iloc[2], 3)), fontsize=8, transform=ax.transAxes, color='red')  # print R2
        ax.text(0.65, 0.48, "lag=" +str(round(stats.iloc[4], 1)), fontsize=8, transform=ax.transAxes, color='red')  # print R2
        # Printing title above plot
    if title != "":  
        ax.set_title(title, fontsize=14, weight="bold")

def detectJump(x,y):
    '''
    Detects if a jump in OD has happened in the slope by taking the derivate.
    The derivate at a jump is significantly different from the average of the 
    last hour of derivates. A jump is thus defined as true if the derivate for 
    the coming point is an outlier compared to the past hour.
    
    Critera for significantly different: >5x st.dev. of average. 
    
    x, y    - Pandas series

    Returns True for positive detection, false otherwise
    '''
    # Convert to DataFrame to remove points with OD below 0.1 (LN = -2.3), where noice is high
    df = pd.concat([x,y], axis=1)
    df = df[df.iloc[:, 1] > -2.3] 
    
    if df.empty:
        return False
    
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]

    # Compare derivates, where the first 7 points represent previous hour
    dt = [] 
    i = 1
    j = 0
    while i+1 < len(x):
        if j < 8: 
            dt.append( (y.iloc[i+1]-y.iloc[i-1])/(x.iloc[i+1]-x.iloc[i-1]) )
        else:
        
            ave = sum(dt[-7:])/len(dt[-7:])
            stdev = np.std(dt[-7:])
            
            newDer1 = (y.iloc[i+1]-y.iloc[i-1])/(x.iloc[i+1]-x.iloc[i-1])
            
            if (newDer1 > (ave + 5 * stdev) ) or ( newDer1 < (ave - 5 * stdev) ):
                return True
            else:
                dt.append(newDer1)
        j += 1
        i += 2
    return False

def getSlopePoints(x, y, w, R2crit):
    '''
    calculates a moving slope in each point of a window specified. Saves the slope if the R2 fulfills 
    the R2 criteria and no jump is detected in the slope. A jump could come from 
    non-biological distrubances in the reading. The slope is saved as an object
    of class SlopePoint together with statistics. Continues to calculate a moving slope.
    x, y    - Pandas series
    w       - integer, minimal window to consider for the slope, specified in span of indeces
    R2crit  - criteria for minimal acceptable R2
    
    Returns list of objects of class SlopePoint
    '''

    slopes = [] 
    i = 0 

    # for jump detection, we need to check one extra point after the window, thus +1 below
    while i < (len(x) - w+1):
        j=i+w

        stats = linregress(x[i:j],y[i:j])
        
        if (stats.rvalue > R2crit):
            if not detectJump(x[i:j+1],y[i:j+1]):
                slopes.append(SlopePoint(stats.slope,i,j,stats.rvalue))
        i += 1
    return slopes

def extendSlope(x,y,maxSlope, R2crit, muCrit):    
    '''
    Extends a slope if it fulfills the criteria: R2 larger than R2crit and mu 
    larger than muCrit of the maximum slope provided in maxSlope. 
    
    x, y       - Pandas series
    maxSlope   - object of SlopePoint containing maximum slope found and statistics
    w          - integer, minimal window to consider for the slope, specified in span of indeces
    R2crit     - criteria for minimal acceptable R2
    muCrit     - criteria for acceptable decrease (0-1) in maximal growth rate
    Returns the extended slope
    '''
    
    extSP = maxSlope
    i = maxSlope.start 
    j = maxSlope.end
    
    upDone = False # boolean for the extention upwards (forward)
    downDone = False # boolean for the extention downwards (backward)
    up = True #boolean for the direction. True: up, False; down

    # calculates new slope and compares to max slope, until extention
    # in both up and down direction is done
    while (not upDone) or (not downDone):
        newSlope = linregress(x[i:j],y[i:j])

        if ((newSlope.slope >= (muCrit * maxSlope.slope)) and (newSlope.rvalue > R2crit)):
            if (not detectJump(x[i:j+1],y[i:j+1])):
                extSP = SlopePoint(newSlope.slope,i,j,newSlope.rvalue)
            else:
                if up:
                    upDone = True
                else:
                    downDone = True
        else:
            if up: #done with extention upwards, go back one step
                upDone = True
                j -= 1
            else: #done with extention downwards, go back one step
                downDone = True
                i += 1
        
        # check 2 points ahead that coming increase/decrease is within boundaries 
        # to ensure jump detection will be possible
        if (j+2 > len(x)):
            upDone = True
        if (i-1 < 0):
            downDone = True        
            
        # switch direction
        up = not up        
        
        # if one direction is done, force the other
        if upDone:
            up = False
        elif downDone:
            up = True

        if up:
            j += 1
        else:
            i -= 1        
    return extSP

def getMaxSlopePoint(x,y, w, R2crit, muCrit):
    '''
    Finds the maximum slope of a window (w) specified and expands it if it 
    fulfills the criteria: R2 larger than R2crit and mu larger than muCrit of 
    the maximum slope found in the first place. Saves the slope as object of class
    SlopePoint, which can carry statistics.
    
    x,y     - Pandas series
    w       - integer, minimal window to consider for the slope, specified in span of indeces
    R2crit  - criteria for minimal acceptable R2
    muCrit  - criteria for acceptable decrease (0-1) in maximal growth rate
    
    Returns extended maxium slope in the class SlopePoint
    '''
    
    #Get list of a moving slope with specified window fulfilling R2 criteria
    slopePoints = getSlopePoints(x,y, w, R2crit) 
    
    # find max point
    maxSP = SlopePoint(0,0,0,0)
    match = False
    
    for sP in slopePoints:
        if (sP.slope > maxSP.slope):
            maxSP = sP
            match = True
    
    # if a slope was found, extend it
    if not match: 
        return maxSP
    else:
        extMaxSP = extendSlope(x, y, maxSP, R2crit, muCrit)
        return extMaxSP

def getStrainPerWell(strainTable):
    '''
    Converts a 8x12 table (2D list), which indicates the position of all strains, 
    into a Pandas DataFrame that lists each well and its strain
    
    strainTable   - a 2D list with 8 x 12 elements containin the strain names at their designated position
    
    Returns Panda DataFrame (1x96 elements) with wells as the columns and one row 
    with the allocated strains
    '''
    strains = []
    wells = []
    rows = ["A","B","C","D","E","F","G","H"]
    
    for i in range(len(rows)):
        for j in range(12):
            wells.append(rows[i]+str(j+1))
            strains.append(strainTable[i][j])
            
    return pd.DataFrame([strains], columns=wells)

def appendToStats(stats, df, text, calculate):
    '''
    Appends the data in df to the variable stats, including text and a calculation
    of average, standard deviation and CV
    stats      - complied list of statistics per strain
    df         - Pandas DataFrame with data to append
    text       - string to append 
    calculate  - Boolean, true if average, SD and CV are to be calculated
    
    Returns the stats variable inputed with appended data
    '''
    stats.append(text)
    for column in df:
        stats.append(df[column].iloc[0])
    if calculate:
        stats.append("Average")
        stats.append(df[df!=0].mean(axis=1).iloc[0])
        stats.append("St dev")
        stats.append(df[df!=0].std(axis=1).iloc[0])
        stats.append("CV (%)")
        # only possible if ave and SD is not zero
        try:
            stats.append(100*(df[df!=0].std(axis=1).iloc[0]/df[df!=0].mean(axis=1).iloc[0]))
        except:
            stats.append("")
    return stats

def main():
    print("Welcome to Plate Reader Analytics!")
    print("After the script is done, start by checking for contamination in blanks. If any is observed, remove them from the list of blanks in the code and rerun the script.")
    print("Extracting data...")
    startTimer = timer()
    
    # Import data
    data_filePath = "OD_data.csv"
    strain_index_filePath = "strain_index.csv"
    media_index_filePath = "media_index.csv"
    strains_filePath = "strain_names.csv"

    if not os.path.exists(data_filePath):
        print("Error: data csv file does not exist")
        sys.exit()
    if not os.path.exists(strain_index_filePath):
        print("Error: 'strain_index_csv' file does not exist")
        sys.exit()
    if not os.path.exists(media_index_filePath):
        print("Error: 'media_index.csv' file does not exist")
        sys.exit()
    if not os.path.exists(strains_filePath):
        print("Error: 'strain_names.csv' file does not exist")
        sys.exit()
    
    # Create directories Analysis, OD and LN-OD, if they don't exist
    if not os.path.exists("Analysis"):
        os.mkdir("Analysis")
    if not os.path.exists("Analysis/OD_plots"):
        os.mkdir("Analysis/OD_plots")
    if not os.path.exists("Analysis/LN-OD_plots"):
        os.mkdir("Analysis/LN-OD_plots")        
        
    ###### Extract media and strain posistions, process it into lists ######
    
    strainsPosTable = pd.read_csv(strain_index_filePath, header=None, engine="python").values.tolist()
    strainPerWell = getStrainPerWell(strainsPosTable)

    mediaPosTable = pd.read_csv(media_index_filePath, header=None).values.tolist()
    strains = pd.read_csv(strains_filePath, header=None).values.tolist()

    # Strains comes as a 2D list from csv --> reformat to 1D list
    strainsReformatted = []
    for s in strains:
        strainsReformatted.append(s[0])
    strains=strainsReformatted
    
    # Media types (here based on carbon sources) which the program will look for
    # and can be set by the user
    media = ['cb','glu','fru']
    
    # Check for blank in strains list
    if "B" not in strains:
        print("Error: missing blank (B) in list of strains. Add and write upper case: B")
        sys.exit()
    
    # make 2D list of indices for each strain/media. The order is according to
    # the strains and media lists above
    strainIndices = getIndicesFromTable(strains, strainsPosTable) 
    mediaIndices = getIndicesFromTable(media, mediaPosTable) 

    if len(strainIndices)==0:
        print("Error: no strains found in strain layout table")
        sys.exit()

    if len(mediaIndices)==0:
        print("Error: no media found in media layout table")
        sys.exit()
    
    ###### Extract OD data and process it for calculations ######
    
    data = pd.read_csv(data_filePath)
    
    #Get time in hours. Save time stamp hh:mm:ss as 3 columns
    time = data["Time"].str.split(":", expand = True)

    time[0]=pd.to_numeric(time[0], errors='ignore') #hours
    time[1]=pd.to_numeric(time[1], errors='ignore') #minutes
    time[2]=pd.to_numeric(time[2], errors='ignore') #seconds

    data["time_h"]=time[0]+time[1].div(60)+time[2].div(3600)
    data=data.drop(columns=["Time","T"])
    
    data.apply(pd.to_numeric, errors='coerce')
    data = data[data.time_h != 0] #removes NaN values at the end

    # Check correct headings
    letters = ["A","B","C","D","E","F","G","H"]
    headings = data.drop(columns = ["time_h"]).columns.values.tolist()
    for i in range(8):
        for j in range(12):
            if headings[i*12+j] != (letters[i] + str(j+1)):
                print("Error! Incorrect formatting of data file. Check headings row.")
                sys.exit()
    
    # find index width which corresponds to specified windowHours
    windowHours = 5
    window = bisect_right(data["time_h"].values, windowHours)
    
    # Subtract the first 1 h due to heating of plate
    data = data[data["time_h"]>1]
    data.reset_index(inplace = True, drop = True) # resets index numbering to start from 0

    # Extract blanks
    blankPos = strainIndices[strains.index("B")]
    dataBlanks = data[blankPos]
    dataBlanks.insert(0,"time_h", data["time_h"])
    
    '''
    # Uncomment this section if blanks need to removed. 
    blanksToRemove = ['A1']    
    dataBlanks = dataBlanks.drop(blanksToRemove, axis=1)
    for i in blanksToRemove:
        blankPos.remove(i)
    '''
    
    #substract a mean of the blanks and then remove them from data
    dataSubBlank = data.subtract(data[blankPos].mean(axis=1), axis='rows')
    dataSubBlank["time_h"]=data["time_h"] 
    data = dataSubBlank.drop(blankPos, axis=1) 

    headings = data.drop(columns = ["time_h"]).columns.values.tolist()

    endTimer = timer()
    print("Took "+str(round(endTimer-startTimer, 0))+" s.")
    print("Calculating growth rates...")
    startTimer = timer()
        
    ######################### Calculation of growth rate #########################    
    dataLN = data.apply(np.log)

    # The quantification limit is LN(OD) = -4.66
    dataLN[dataLN < (-4.66)] = np.NaN
    # The linearity extends up to LN(OD) = -0.5
    dataLN[dataLN > (-0.5)] = np.NaN
    dataLN["time_h"]=data["time_h"]
    
    #criteria for the slope: R2 >= 0.995, mu >= 95% of max slope
    R2crit = 0.995
    muCrit = 0.95
    
    # Find highest slope
    maxSlopePoints = []
    print("Done with well: ")
    for column in headings:
        maxSlopePoints.append(getMaxSlopePoint(dataLN["time_h"], dataLN[column], window, R2crit, muCrit))
        print(column),
    print("")
    
    # Extract the values
    growthRates = []
    r_values = []
    range_mu = []
    start_pos = [] #of exponential phase, indicated as index
    end_pos = [] #of exponential phase, indicated as index
    diff_time = data["time_h"].iloc[2]-data["time_h"].iloc[1] 
    
    for i in range(len(maxSlopePoints)):
        growthRates.append(maxSlopePoints[i].slope)
        r_values.append(maxSlopePoints[i].r_value)
        range_mu.append((maxSlopePoints[i].end-maxSlopePoints[i].start)*diff_time)
        start_pos.append(maxSlopePoints[i].start)
        end_pos.append(maxSlopePoints[i].end)

    # convert all lists to DataFrames 
    growthRates = pd.DataFrame([growthRates], columns = headings)
    r_values = pd.DataFrame([r_values], columns = headings)
    range_mu = pd.DataFrame([range_mu], columns = headings)
    start_pos = pd.DataFrame([start_pos], columns = headings)
    end_pos = pd.DataFrame([end_pos], columns = headings)
    
    endTimer = timer()
    print("Took " + str(round(endTimer - startTimer, 0)) + " s.")
    print("Calculating remaining parameters and statistics...")
    startTimer = timer()
    
    ######################### Get max OD #########################
    data_without_time = data.drop(columns = ["time_h"])
    maxOD = data_without_time.max()
    maxOD = pd.DataFrame(maxOD).transpose()
    
    ######################### Get lag time #########################
    # Lag time defined as start of exponential phase, thus incl. acceleration phase
    lag_times = []
    for column in headings:
        lag_times.append(data['time_h'].iloc[start_pos[column]].iloc[0]) 
    lag_times = pd.DataFrame([lag_times], columns = headings)
    
    endTimer = timer()
    print("Took "+str(round(endTimer-startTimer, 0))+" s.")
    print("Compiling statistics into csv files...")
    startTimer = timer()
    
    ######################### Generate a statistics table #########################
    # logging of calculated data into csv files
    # compile into one table, reporting stastistics per well
    statsTable = pd.DataFrame(columns = headings)
    statsTable = statsTable.append(strainPerWell, ignore_index = True)
    statsTable = statsTable.append(growthRates, ignore_index = True)
    statsTable = statsTable.append(r_values, ignore_index = True)
    statsTable = statsTable.append(maxOD, ignore_index = True)
    statsTable = statsTable.append(lag_times, ignore_index = True)
    statsTable = statsTable.append(start_pos, ignore_index = True)
    statsTable = statsTable.append(end_pos, ignore_index = True)
    statsTable = statsTable.append(range_mu, ignore_index = True)
    statsTable = statsTable.transpose()
    statsTable.to_csv('Analysis/statistics_raw.csv', header = ['Strain', 'Growth rate (1/h)', 'R2', 
                                                      'Max OD', 'Lag time (h)', 'Start exp point',
                                                      'End exp point', 'Range (h)',
                                                      ]) # save as csv
    
    #create another table with compiled statistics per strain and media
    statsCompiled=[]
    for i in range(len(media)):
        for j in range(len(strains)): 
            s=strains[j]
            pos = intersection(mediaIndices[i],strainIndices[j])
            stats = []
            if len(pos) > 0 and not strains[j]=="B":
                stats.append(strains[j])
                stats.append(media[i])
                for p in pos:
                    stats.append(p)
                stats = appendToStats(stats, growthRates[pos], "Growth rate (1/h)", True)
                stats = appendToStats(stats, r_values[pos], "R2", False)
                stats = appendToStats(stats, lag_times[pos], "Lag time (h)", True)
                stats = appendToStats(stats, maxOD[pos], "Max OD", True)
                stats = appendToStats(stats, range_mu[pos],"Range (h)", False)
                stats = appendToStats(stats, start_pos[pos],"Start exp point", False)
                stats = appendToStats(stats, end_pos[pos],"End exp point", False)
                statsCompiled.append(stats)
    statsCompiled = pd.DataFrame(statsCompiled).transpose() 
    statsCompiled.to_csv('Analysis/statistics_compiled.csv', header=False, index=False) # save as csv
    
    endTimer = timer()
    print("Took " + str(round(endTimer - startTimer, 0)) + " s.")
    print("Plotting...")
    startTimer = timer()

    #########################  Plotting #########################
    # Plot both OD data and LN(OD) data. 
    dataLN = data.apply(np.log)
    dataLN["time_h"]=data["time_h"]
    dataLN.to_csv('Analysis/lnDataForPlotting.csv') # save as csv
    
    # remove values below detection limit (LN(OD) = -6.0)
    dataLN[dataLN<-6.0] = np.NaN
    
    print("Done with: ")
    statsTable = statsTable.transpose() 
    for i in range(len(strains)):
        pos = []
        for j in range(len(mediaIndices)): 
            pos.append(intersection(mediaIndices[j],strainIndices[i]))
        
        print (strains[i]),
        if(strains[i] == "B"):
            plotStrain(strains[i], dataBlanks, statsTable, pos, media, "Analysis/OD_plots/Growth curves for blanks.pdf")
        else:
            plotStrain(strains[i], data, statsTable, pos, media, "Analysis/OD_plots/Growth curves for "+strains[i]+".pdf", 
                          draw_slopePoints = True, draw_lag_time = True, print_stats = True)
            plotStrain(strains[i], dataLN, statsTable, pos, media, "Analysis/LN-OD_plots/Growth curves for "+strains[i]+" - LN.pdf", 
                          plotLN = True, draw_slopePoints = True, draw_lag_time = True, print_stats = True)
    
    endTimer = timer()
    print("")
    print("Took " + str(round(endTimer - startTimer, 0)) + " s.")
    print("Program finished.")

if __name__== "__main__":
  main()