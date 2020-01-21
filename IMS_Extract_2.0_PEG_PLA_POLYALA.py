#!/usr/bin/python

#------------
#
#IMS Extract v2.0
#
#Written by qduez - November 2019
#Updated by shoyas - January 2020
#Licensed under GNU General Public License v3.0; See license file.
#------------
#Requires Python 3 or later (tested with 3.7.5) and TWIMExtract (https://sites.lsa.umich.edu/ruotolo/software/twim-extract/)

#PLEASE CITE THESE PAPERS WHEN USING THIS SCRIPT: https://doi.org/10.1021/acs.analchem.7b00112 (TWIMExtract), https://doi.org/10.1007/s13361-017-1762-4 (PEG, PLA calibration), https://doi.org/10.1021/ac3014498 (PolyALA calibration)

#This script is intended to be used for Traveling Wave Ion Mobility Spectrometry (TWIMS) calibration, with synthetic polymers as calibrants.

#Polymers taken into account: Poly(ethylene glycol) (PEG), Polylactide (PLA) and Polyalanine (POLYALA).

#Polymer sizes taken into account: Commercial PEG (Mn = 600, 1000, 2000 and 3350), PLA (Mn = 4000, 5500), commercial POLYALA (Mn ~ 2000)

#For more info on the polymers used for this calibration : https://doi.org/10.1007/s13361-017-1762-4 (PEG, PLA - Open access) and https://doi.org/10.1021/ac3014498 (PolyALA)

#The aim of this script is to extract the ATDs of polymer ions with givens masses (or DPs) for charge states between 1+ and 3+.

#This script will extract the ATDs from Waters .raw files for the mass ranges corresponding to PEG, PLA and PolyALA ions across all spectra recorded by the user
#You need this script, GaussFit and your MS spectra in the same folder

#The extraction will be performed automatically using TWIMSExtract developped by the Ruotolo Lab (10.1021/acs.analchem.7b00112)

#Each ion familly is then to be processed by GaussFit which fits each ATD by a Gaussian function and yields the maximum of each fitted curve. The fitted data can then be used to calibrate TWIMS data.

#!!!! THE WINDOWS PATH OF INSTALLATION OF TWIMSExtract SHOULD BE ADDED HERE BETWEEN THE "":

twexpath = (r"C:\TWIMExtract\jars")

import os
import subprocess
import re
import csv
import glob
import pandas as pd
import fnmatch
import shutil

#Precise masses for elements in PEG, PLA and PolyALA (As defined in MassLynx)

C=12.0000
H=1.0078
O=15.9949
Na=22.9898
N=14.0031

PEG=(2*C+4*H+O)
PLA=(6*C+8*H+4*O)
PolyALA=(3*C+5*H+O+N)

pwd=os.getcwd()

charge=[1,2,3]

# Definition of functions to: 
# 1) create TWIMExtract "range files" containing mass range for PEG, PLA and PolyALA. These files will further be used when calling TWIMextract

# 2) make use of fuction #1 for the data extraction with specific conditions depending on the polymer size and charge

# 3) call TWIMExtract recursively on all "range files" created by #2

# 4) get the number of lines

# 5) recreate an arrival time scale based on pusher frequencies recorded in the '_extern.inf' file 

def  make_mass_filesPEG(num1,num2,mass1,mass2,Name):
# Function that will create the formatted range files for TWIMSExtract (for PEG)
# Here, num1 is the charge and num2 is the DP and Name corresponds to the polymer name contained in the '_HEADER.txt' file
    os.chdir("Input_" + Name + "_" + str(num1) + "+")
    f = open("PEG_" + str(num2).zfill(2) + "_z" + str(num1) + ".txt","w+")
    f.write("MZ_start_(m/z): " + str(mass1) + "\n")
    f.write("MZ_end_(m/z): " + str(mass2) + "\n")
    f.write("RT_start_(minutes): 0" + "\n")
    f.write("RT_end_(minutes): 100" + "\n")
    f.write("DT_start_(bins): 1" + "\n")
    f.write("DT_end_(bins): 200")
    f.close()
    os.chdir('..')
def  make_mass_filesPOLYALA(num1,num2,mass1,mass2,Name):
# Function that will create the formatted range files for TWIMSExtract (for POLYALA)
# Here, num1 is the charge and num2 is the DP and Name corresponds to the polymer name contained in the '_HEADER.txt' file
    os.chdir("Input_PolyALA_" + str(num1) + "+")
    f = open("PolyALA_" + str(num2).zfill(2) + "_z" + str(num1) + ".txt","w+")
    f.write("MZ_start_(m/z): " + str(mass1) + "\n")
    f.write("MZ_end_(m/z): " + str(mass2) + "\n")
    f.write("RT_start_(minutes): 0" + "\n")
    f.write("RT_end_(minutes): 100" + "\n")
    f.write("DT_start_(bins): 1" + "\n")
    f.write("DT_end_(bins): 200")
    f.close()
    os.chdir('..')
def  make_mass_filesPLA(num1,num2,mass1,mass2,Name):
# Function that will create the formatted range files for TWIMSExtract (for PLA)
# Here, num1 is the charge and num2 is the DP and Name corresponds to the polymer name contained in the '_HEADER.txt' file
    os.chdir("Input_" + Name + "_" + str(num1) + "+")
    f = open("PLA_" + str(num2) + "_z" + str(num1) + ".txt","w+")
    f.write("MZ_start_(m/z): " + str(mass1) + "\n")
    f.write("MZ_end_(m/z): " + str(mass2) + "\n")
    f.write("RT_start_(minutes): 0" + "\n")
    f.write("RT_end_(minutes): 100" + "\n")
    f.write("DT_start_(bins): 1" + "\n")
    f.write("DT_end_(bins): 200")
    f.close()
    os.chdir('..')

def cali(Name):
# This function will check the type of polymer and apply different conditions depending on the polymer size and mass
# Naming corresponds to the polymer name contained in the '_HEADER.txt' file
    for z in charge:
        #First we prepare the folder for the range files
        f = open(str(z) + '+_' + Name +'_Masses.txt', 'w+')
        f.write ("DP"+ "\t" + "Mass" + "\n")
        f.close()
    
        try:
                os.mkdir("Input_" + Name + "_" + str(z) + "+")
        except OSError:
                pass
        # Initialize the current working directory (here in the polymer directory)
        pwd2 = os.getcwd()    

        # Depending on the polymer, different strategies are employed
        # Example:
        # 1) For PEG with charge state +1, PEG3350 is not taken into account (Larger polymer will not likely be singly charged)
        # 2) For PEG with charge state +3, PEG600 is not taken into account (Shorter polymer will not likely be triply charged)
        
        #Also, only given DPs will be considered (DP8 to 37 for PEG 1+, DP 13 to 26 and 51 to 82 for PEG 2+). For more info, see https://doi.org/10.1007/s13361-017-1762-4
        if "PEG" in Name:
            print(Name+" is PEG.")
            if z == 1:
                print("Charge state +1:")
                i = 8
                # We then create the range files for the extraction using make_mass_filesPEG
                while i <= 37:
                    mass = ((i*PEG + z*Na + (2*H+O))/z) 
                    f = open(str(z) + '+_' + Name +'_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPEG(z,i,start_mass,end_mass,Name)
                    i += 1
                    
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws = glob.glob("*.raw")
                print("Extracting data for " + raws[0] + "\n")
                # For PEG with charge state +1, PEG3350 is not taken into account: raws list will be empty if PEG3350
                with open (raws[0]+"/_HEADER.TXT","r") as f:
                    word = ["PEG3350","PEG 3350","3350"]
                    file_contents = f.readlines()
                    # Line 11 (10 in Python) corresponds to the comment introduced during the measurement
                    my_line = file_contents[10]
                    my_line = my_line.split()
                    result = any(elemen in word for elemen in my_line)
                    if result:
                        raws = []
                if not raws:
                    print('Calibration file contains PEG3350 and will not be taken into account for charge state +1.')
                else:
                    for f in files:
                        print(f)
                        call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                        
                    for f in files:
                        dplist = re.findall("PEG_(\d+)_z" + str(z) +".txt",f)
                        dp=dplist[0]
                        i = 0
                        os.chdir(pwd2+"/Output")
                        # UPDATE: since we do not merge ATDs anymore, we can just rename the files with the extension '.out' for GaussFit
                        os.rename("DT_" + str(raws[i]) +"_fn-1_#PEG_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PEG_z" + str(z) + "_DP" + str(dp).zfill(2) + ".txt_raw.out")
                        os.chdir("..")
                        
            elif z == 2:
                print("Charge state +2:")
                # Two DP range to avoid the folding region
                i = 13
                while i <= 26:
                    mass=((i*PEG + z*Na + (2*H+O))/z)
                    f=open(str(z) + '+_' + Name + '_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPEG(z,i,start_mass,end_mass,Name)
                    i += 1
            
                i = 51
                while i <= 82:
                    mass = ((i*PEG + z*Na + (2*H+O))/z)
                    f = open(str(z) + '+_' + Name + '_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPEG(z,i,start_mass,end_mass,Name)
                    i += 1
                
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws = glob.glob("*.raw")
                print("Extracting data for " + raws[0] + "\n")
                #For each raw file and for each rule file, we make the extraction by calling call_TWIMExtract
                for f in files:
                    print(f)
                    call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                
                #We now have an extracted file for each raw file, we thus need to sum all the extracted files for each raw file
            
                for f in files:
                    dplist = re.findall("PEG_(\d+)_z" + str(z) +".txt",f)
                    dp = dplist[0]
                    i = 0
                    os.chdir(pwd2+"/Output")
                    os.rename("DT_" + str(raws[i]) +"_fn-1_#PEG_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PEG_z" + str(z) + "_DP" + str(dp) + ".txt_raw.out")
                    os.chdir('..')

            elif z == 3:
                print("Charge state +3:")
                i = 32
                while i <= 56:
                    mass=((i*PEG + z*Na + (2*H+O))/z)
                    f=open(str(z) + '+_' + Name + '_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPEG(z,i,start_mass,end_mass,Name)
                    i += 1
                    
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws =  glob.glob("*.raw")
                print("Extracting data for " + raws[0] + "\n")
                # For PEG with charge state +3, PEG600 is not taken into account: raws list will be empty if PEG600
                with open (raws[0]+"/_HEADER.TXT","r") as f:
                    word = ["PEG 600", "PEG600","600"]
                    file_contents = f.readlines()
                    my_line = file_contents[10]
                    my_line = my_line.split()
                    result = any(elemen in word for elemen in my_line)
                    if result:
                        raws = []
                if not raws:
                    print('Calibration file contains PEG600 and will not be taken into account for charge state +3.')
                else:
                    for f in files:
                        print(f)
                        call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                        
                    for f in files:
                        dplist = re.findall("PEG_(\d+)_z" + str(z) +".txt",f)
                        dp=dplist[0]
                        i = 0
                        os.chdir(pwd2+"/Output")
                        os.rename("DT_" + str(raws[i]) +"_fn-1_#PEG_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PEG_z" + str(z) + "_DP" + str(dp) + ".txt_raw.out")
                        os.chdir('..')
                
        elif "POLYALA" in Name:
            print(Name + " is PolyALA.")
            if z == 1:
                print("Charge state +1:")
                i = 3
                #We then create the range files for the extraction using make_mass_files 
                while i <= 14:
                    mass = ((i*PolyALA + z*H + (2*H+O))/z)
                    f = open(str(z) + '+_' + Name +'_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPOLYALA(z,i,start_mass,end_mass,Name)
                    i += 1
                    
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws = glob.glob("*.raw")
                print("Extracting data for " + raws[0] + "\n")
                for f in files:
                    print(f)
                    call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                        
                for f in files:
                    dplist = re.findall("PolyALA_(\d+)_z" + str(z) +".txt",f)
                    dp=dplist[0]
                    i = 0
                    os.chdir(pwd2+"/Output")
                    os.rename("DT_" + str(raws[i]) +"_fn-1_#PolyALA_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PolyALA_z" + str(z) + "_DP" + str(dp).zfill(2) + ".txt_raw.out")
                    os.chdir('..')

            elif z == 2:
                print("Charge state +2:")
                i = 11
                while i <= 26:
                    mass=((i*PolyALA + z*H + (2*H+O))/z)
                    f=open(str(z) + '+_PolyALA_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPOLYALA(z,i,start_mass,end_mass,Name)
                    i += 1
                    
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws = glob.glob("*.raw")
                #For each raw file and for each rule file, we make the extraction by calling call_TWIMExtract
                for f in files:
                    print(f)
                    call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                #We now have an extracted file for each raw file, we thus need to sum all the extracted files for each raw file
                for f in files:
                    dplist = re.findall("PolyALA_(\d+)_z" + str(z) +".txt",f)
                    dp = dplist[0]
                    i = 0
                    os.chdir(pwd2+"/Output")
                    os.rename("DT_" + str(raws[i]) +"_fn-1_#PolyALA_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PolyALA_z" + str(z) + "_DP" + str(dp) + ".txt_raw.out") 
                    os.chdir('..')
            
            elif z == 3:
                print("Charge state +3:")
                i = 19
                while i <= 33:
                    mass=((i*PolyALA + z*H + (2*H+O))/z)
                    f=open(str(z) + '+_PolyALA_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPOLYALA(z,i,start_mass,end_mass,Name)
                    i += 1
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws = glob.glob("*.raw")
                #For each raw file and for each rule file, we make the extraction by calling call_TWIMExtract
                for f in files:
                    print(f)
                    call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                for f in files:
                    dplist = re.findall("PolyALA_(\d+)_z" + str(z) +".txt",f)
                    dp = dplist[0]
                    i = 0
                    os.chdir(pwd2+"/Output")
                    os.rename("DT_" + str(raws[i]) +"_fn-1_#PolyALA_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PolyALA_z" + str(z) + "_DP" + str(dp) + ".txt_raw.out")
                    os.chdir('..')
                    
        elif "PLA" in Name:
            print(Name +" is PLA.")
            if z == 2:
                print("Charge state +2:")
                i = 22
                #We then create the range files for the extraction using make_mass_files 
                while i <= 37:
                    mass = ((i*PLA + z*Na + (C+4*H+O))/z) 
                    f = open(str(z) + '+_PLA_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPLA(z,i,start_mass,end_mass,Name)
                    i += 1
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws = glob.glob("*.raw")
                #For each raw file and for each rule file, we make the extraction by calling call_TWIMExtract
                for f in files:
                    print(f)
                    call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                #We now have an extracted file for each raw file, we thus need to sum all the extracted files for each raw file
                for f in files:
                    dplist = re.findall("PLA_(\d+)_z" + str(z) +".txt",f)
                    dp = dplist[0]
                    i = 0
                    os.chdir(pwd2+"/Output")
                    os.rename("DT_" + str(raws[i]) +"_fn-1_#PLA_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PLA_z" + str(z) + "_DP" + str(dp) + ".txt_raw.out")
                    os.chdir('..')
            
            elif z == 3:
                print("Charge state +3:")
                i = 36
                while i <= 50:
                    mass=((i*PLA + z*Na + (C+4*H+O))/z)
                    f=open(str(z) + '+_PLA_Masses.txt', 'a+')
                    f.write(str(i)+ "\t" + str(mass)+ "\t" + "\n")
                    f.close()
                    end_mass = mass+((0.1)/z)
                    start_mass = mass-((0.1)/z)
                    make_mass_filesPLA(z,i,start_mass,end_mass,Name)
                    i += 1
                    
                files = os.listdir("./Input_" + Name + "_" + str(z) + "+/")
                raws = glob.glob("*.raw")
                #For each raw file and for each rule file, we make the extraction by calling call_TWIMExtract
                for f in files:
                    print(f)
                    call_TWIMExtract(z,raws[0],f,pwd2,PolyName)
                #We now have an extracted file for each raw file, we thus need to sum all the extracted files for each raw file
                for f in files:
                    dplist = re.findall("PLA_(\d+)_z" + str(z) +".txt",f)
                    dp = dplist[0]
                    i = 0
                    os.chdir(pwd2+"/Output")
                    os.rename("DT_" + str(raws[i]) +"_fn-1_#PLA_"+ str(dp) +"_z" +str(z) + ".txt_raw.txt", "ATD_PLA_z" + str(z) + "_DP" + str(dp) + ".txt_raw.out")
                    os.chdir('..')
                    
def call_TWIMExtract(num1,rw,rng,folder,Name):
    #Function to call TWIMExtract through Java
    #Here, num1 is the charge, rw is the location of the raw file and rng is the name of the range file
    # Name is the polymer name contained in the '_HEADER.txt' file
    subprocess.call(['java', '-jar',str(twexpath) + '\TWIMExtract.jar','-i', str(folder) + "/" + str(rw),'-o', str(folder) +'/Output/','-m','1','-r',str(folder) + '/Input_' + Name + '_' + str(num1) + '+/' + str(rng)])
    os.chdir("Output/")
    #Deletes the first lines of the csv file
    with open('DT_'+str(rw)+'_fn-1_#'+str(rng)+'_raw.csv', 'r') as fin:
        data = fin.read().splitlines(True)
    with open('DT_'+str(rw)+'_fn-1_#'+str(rng)+'_raw.csv', 'w') as fout:
        fout.writelines(data[3:])
    #Converts the file from a "comma separated file" to a tab-separated file
    with open('DT_'+str(rw)+'_fn-1_#'+str(rng)+'_raw.csv','r') as fin:
        with open('DT_'+str(rw)+'_fn-1_#'+str(rng)+'_raw.txt', 'w',newline='') as fout:
            reader = pd.read_csv(fin)
            reader.drop(reader.columns[0], axis = 1, inplace = True)
            reader.insert(0,column='A',value = driftlist)
            reader.to_csv(fout, sep='\t',header = False, index = False)
    os.chdir('..')
    

def file_len(fname):
    #Function to get the line count of a file ; From https://stackoverflow.com/questions/845058/how-to-get-line-count-of-a-large-file-cheaply-in-python
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def corr_drif (rwf):
# Will re-create an arrival time scale based on the '_extern.inf' file contained in the first spectrum in the folder

# We do not use the arrival times provided by TWIMExtract since they are susceptible to have a slight shift (max 0.1 ms) for Synapt G2-S or G2-Si.
# This will rereate a list of real arrival times!
# It will read the pusher frequency in the '_extern.inf' file, add the offset (0.25 µs) and multiply by the bin range.

    with open(rwf+"/_extern.inf", 'rb') as f:
         line = f.readline()
         cnt = 1
         while line:
             if b'ADC Pusher Frequency' in line:
                 line = line.decode('latin-1')
                 Pushf = line.split()[4]
                 PushF = float(Pushf) + 0.25
                 PushF = PushF/1000
                 bstart = 1
                 with open("driftimes", 'w') as test:
                     while bstart <= 199:
                         time = round(bstart * PushF,3)
                         test.write(str(time)+'\n')
                         bstart += 1
             line = f.readline()
             cnt += 1

# List raw files
RawList = glob.glob("*.raw")

corr_drif(RawList[0])
drift = open("driftimes")
driftlist = drift.read().splitlines()

SpecList = []
# Start the actual program
for file in RawList:
    print(file)
    f = open(file+"/_HEADER.txt", 'r')
    i = f.readlines()
    my_line = i[10]
    PolyName = "P"+str(my_line.split("P",1)[1])
    PolyName = PolyName.replace(" ","")
    PolyName = PolyName.replace("\n","")
    PolyName = PolyName.upper()
    SpecList.append(PolyName)
    print("Name found in _HEADER.txt file: " + PolyName)
    try:
        os.mkdir(PolyName)
        print('Directory called "' + PolyName + '" created.')
        shutil.copytree(pwd+"/"+file, pwd+"/"+PolyName+"/"+file)
        print('MS spectrum ' + file + ' copied in ' + PolyName + '.')
    except OSError:
        pass
    os.chdir(pwd+"/"+PolyName+"/")
    print("Go to " + pwd + "/"+ PolyName + "/ ...")
    try:
        os.mkdir("Output")
    except OSError:
        pass
    cali(PolyName)
    
    #Calls GaussFit
    print("Starting GaussFit...")
    subprocess.Popen('python ../GaussFit_0.4i_3.7.5.py > Output.log', shell=True).wait()
    print("Over!")
    shutil.copyfile(str(pwd)+"\\"+PolyName+"\\Output.log", pwd+"\\"+PolyName+".log")
    print("GaussFit output copied as " + PolyName+".log" + " in " + pwd)
    os.chdir("..")
    
print("Polymers used for calibration: " + str(SpecList))
