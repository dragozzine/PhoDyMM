# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:13:50 2021

@author: Daniel
"""

import numpy as np
import os
import shutil
import glob 
import itertools
import csv
import matplotlib.pyplot as plt
import math

RUNSPATH = "C:/Users/Daniel/Documents/Graduate_Research/PhoDyMM-master/testedSystems/Test5_modelAnalysis/"
WATER = 'Purewater_MvR.dat'
HE = '1per_HHe_MvR.dat'
EARTH = 'Earthlike_MvR.dat'
PLANETCOLS = 9 #number of columns of data for each planet
MASSCOL = 7 #mass column within each set of nine columns
RADIUSCOL = 8 #radius column within each set of nine columns

RSUN = 696340000
MJUP = 1.898 * 10**27
RJUP = 69911000
MEARTH = 5.972 * 10**24
REARTH = 6371000
maxMass = 20
maxRadius = 10

plt.rcParams.update({'font.size': 45})

#Finds the names of all the systems in the runs file
def getSystems():
    
    systems = [ f.name for f in os.scandir(RUNSPATH) if f.is_dir() ]
    
    return systems

#Using the Parameters File Header, count how many planets there are 
def numPlanets(header):
    return header.count('Planet')

#Reads in the last 100 lines of the dqa file for each system, saves mass and radius for each planet
    #returns mass/massJupiter,radius in meters all planets
def parseData(systems):
    
    mass = []
    radius = []
    lineToRead = 100
    numSystems = len(systems)
    #numSystems = 1
    
    for i in range(numSystems):
        fileName = RUNSPATH + systems[i] + '/analysis_dir/dqa_allparam.csv' 
        #print(fileName)
        try:
            file = open(fileName)
            header = next(file)
            nPlanets = numPlanets(header)
            rstarIn = PLANETCOLS * nPlanets + 2
            
            for line in (file.readlines()[-lineToRead:]):
                line = line.split(',')
                
                for j in range(nPlanets):
                    mindex = 1+ MASSCOL + PLANETCOLS * j
                    rindex = 1+ RADIUSCOL + PLANETCOLS * j
                    mass.append(float(line[mindex])) #measured in m/mjup
                    radius.append(float(line[rindex]) * float(line[rstarIn]) *RSUN) #measured in meters
                    
            
            
        except OSError as error:
            print(error)
               
        
    return mass,radius

#reads in mass radius distributions given in data files   
def readDistributions(fileName):
    mass = []
    radius = []
    file = open(fileName)
    
    for line in file:
        line = line.split()
        mass.append(float(line[0]))
        radius.append(float(line[1]))
    
    file.close()
    return mass,radius

#given a set density and an array of mass values, returns the radius values
def calcIsoRadius(mass,rho):
    radius = []
    
    for m in mass:
        mSI = m * MEARTH
        rSI = math.pow(3*mSI/(4*rho*math.pi),1/3)
        rE = rSI/REARTH
        radius.append(rE)
    
    return radius
        

def makeMassRadiusDis():
    systems = getSystems()
    mass,radius = parseData(systems)
    
    mass = np.array(mass)
    radius = np.array(radius)
    
    mass = mass*MJUP/MEARTH #Convert to Earth units
    radius = radius/REARTH
    
    waterBlack = readDistributions(WATER) #get specific distributions
    hHeEarthRed = readDistributions(HE)
    earthGreen = readDistributions(EARTH)
    
    isoMass = np.linspace(0,maxMass,1000) #find isodensity values
    r1Iso = calcIsoRadius(isoMass,1000)
    r30Iso = calcIsoRadius(isoMass,30000)

    plt.scatter(mass,radius,alpha=.15)
    plt.xlabel('Mass/Mass Earth')
    plt.ylabel('Radius/Radius Earth')
    plt.ylim([0,maxRadius])
    plt.xlim([0,maxMass]) 
    plt.plot(waterBlack[0],waterBlack[1],'k',label='Water World',linewidth = 7)
    plt.plot(hHeEarthRed[0],hHeEarthRed[1],'r', label = '1% H/He',linewidth = 7)
    plt.plot(earthGreen[0],earthGreen[1],'g',label = 'Earth-like',linewidth = 7)
    plt.plot(isoMass,r1Iso,'-b',label = 'rho = 1',linewidth = 7)
    plt.plot(isoMass,r30Iso,'-m', label = 'rho = 30',linewidth = 7)
    plt.xlim([0,10])
    plt.ylim([0,6])
    ax = plt.gca() #you first need to get the axis handle
    ax.set_aspect(1.4) #sets the height to width ratio to 1.5. 
    plt.show()