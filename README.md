# Plate reader analytics

Table of Contents
-----------------
1. About the script
2. Getting started (prerequisits, installation)
3. Usage
4. Contact

About the script
----------------
v7.1 - 2020-12-08  
This script was developed by Johannes Yayo.

The aim with this script is to analyze OD data generated in a plate reader. The script calculates growth parameters (growth rate, lag time, max OD) including statistics as well as plots OD and LN(OD). It is specifically adapted to data generated from the Biotek Epoch 2 plate reader. It is coded in Python 2.7 (from the Anaconda distribution). 

The script is written for an experiment where several strains are tested in different media. Blanks are expected and will be substracted from the OD data at each time point. If a blank was contaminated it can be excluded from the list of blanks in the code as indicated. It is recommended to inspect each plot for misfit data and outliers. 

The growth rate is defined as the slope to a linear fit to the natural logaritm of OD. By default, the fit should have a R2 > 0.995 over a window of 5 hours, representing the exponential phase.

For the plate reader specified above, dilution series were used to find its detection limit, quantification limit and linearity. It was based on Clostridium thermocellum cells harvested at mid-exponential phase. Based on those, following limits were found (expressed in plate reader OD):  
* detection limit = 0.0025  
* quantification limit = 0.01  
* linear up to 0.6 

Getting started
---------------

Prerequisits:
* Python 2.7
* OD data from a 96-well plate
* Blanks are included

Installation
* Download the script and place it in the same folder as the OD data

Usage
-----
Preparation:
1.	Export OD data from the plate reader as an excel file. Each column must represent each well, whereas the rows represent time points. Include all 96 wells and temperature. Place the time-stamp column first, then temperature and then all wells in order (A1, ..., A12, B1, ..., H12). Rename the temperature column to “T” and the time column "Time". Leave the time in the format HH:MM:SS. Convert the excel file into a csv file and save as 'OD_data.csv'. 
2.	Make three additional csv-files:  
i. list of strain names tested including blank ('B'). Save as 'strain_names.csv'.
ii. strain positions in the plate layout. Consists of a plain 8 x 12 (rows x columns) table representing the 96-well plate (without header or special first column). Each position in the table/layout is indicated with a strain name or blank ('B'). Save file as 'strain_index.csv'.  
iii. media type in the plate layout. Consists of a plain 8 x 12 (rows x columns) table representing the 96-well plate (without header or special first column). Each position in the table/layout is indicated with a media type. In its current state, the script will look for media named cb, glu, or fru, and at least one need to be present. For other media, adjust in the script. Save file as 'media_index.csv'.  
3.	Place all files in the same directory as the scipt

Run the plate_reader_analytics.py file either in command line or in a Python environment (e.g. Spyder from Anaconda distribution).

The program generates a folder called Analysis. This folder contains OD plots, LN(OD) plots and data-files. The LN(OD) data is also exported and can be used for manual curation. Note that the first hour and OD values below the detection limit is removed.

Post-run checks:
1.	Contamination in blanks. Open the file “Growth curves for blanks.pdf” in the folder OD_plots. If contamination is seen, remove that well from the blanks. To do this, open the script for editing (e.g. in Spyder or Notepad). In the main function, find the block seen below and replace A1 with the blank to remove. Then remove the ''' marks before and after the block.
>'''   
>\# Uncomment this section if blanks need to removed.  
>blanksToRemove = ['A1']  
>dataBlanks = dataBlanks.drop(blanksToRemove, axis=1)  
>for i in blanksToRemove:  
>    blankPos.remove(i)  
>'''  
2. Rerun the script if blanks were removed.
3. (Recommended) Open the folder “LN-OD_plots”. Go through each growth curve and check that the points included for calculating the growth rate (red points) cover a relevant part of the exponential phase. Note any erroneous jumps or other non-biological trends in OD that can influence the calculation. 
4. Further analysis:
The file "statistics_raw.csv" outputs a table with parameters and statistics collected per well
The file “statistics_compiled.csv” outputs a table with growth parameters and statistics per strain and media type.

Contact
-------
Johannes Yayo
johannes.yayo@gmail.se
