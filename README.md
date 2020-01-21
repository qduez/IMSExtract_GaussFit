# IMSExtract & GaussFit
Python script that automatically extracts and fits Traveling Wave Ion Mobility Mass Spectrometry data for calibration purposes
Authors: Sébastien Hoyas & Quentin Duez (Organic Synthesis and Mass Spectrometry Laboratory - UMONS https://s2mos.umons.ac.be/)

Feel free to contact us : sebastien.hoyas@umons.ac.be / quentin.duez@umons.ac.be

These scripts are licensed under GNU General Public License v3.0, see License file.

When using these scripts for your research, please cite: https://doi.org/10.1021/acs.analchem.7b00112 (TWIMExtract), https://doi.org/10.1007/s13361-017-1762-4 (PEG, PLA calibration), https://doi.org/10.1021/ac3014498 (PolyALA calibration)

1. Preamble

These scripts should be used in the frame of Traveling Wave Ion Mobility Spectrometry calibration with polymer ions as calibrants. 
Recently, synthetic polymers were introduced as new TWIMS calibrants. Due to the large amount of calibrating ions, the manual extraction of IMS data can be cumbersome. 
The scope of these scripts is to automatically extract arrival time distributions (ATD) of polymers ions covered by the polymer calibration and the polyalanine calibration. 
This script only works on .raw files recorded with Waters instruments.

The extracted ATDs are fitted with Gaussian functions by the script GaussFit. Maxima are then reported and can be used to calibrate TWIMS experiments. 
This script is made to use recursively TWIMExtract, developed by the Ruotolo Lab.

The purpose of these scripts is to make the extraction of data for the calibration easier, not to replace manual extraction. 
Fitting errors can occur depending on the way you record your spectra or on your spectra quality.

2. Requirements

- Python 3.7.5
Installation (Windows) - For Linux see GaussFit script:
Install Python on Windows: https://www.python.org/ftp/python/3.7.5/python-3.7.5rc1-amd64.exe
Add Python to PATH variable
Open Windows command prompt
Install PIP: https://www.liquidweb.com/kb/install-pip-windows/
Install modules:
pip install lmfit
pip install pandas
pip install matplotlib

- TWIMExtract (https://sites.lsa.umich.edu/ruotolo/software/twim-extract/)

!!! Include TWIMExtract installation path in the beginning of the IMSExtract script !!!
Note that Java is required for TWIMExtract to work.


3. Usage

- Record your calibration spectra with a Waters instrument and specify the name of the polymer (PEG, PLA or POLYALA) with the approximate Mn in the spectrum header.
!!! This is mandatory as the script applies different conditions depending on the polymer type and size !!!

- Create a new folder for the execution of the scripts.
!!! Don't execute the script twice in the same folder without deleting the previously generated data !!!

- Copy-paste the scripts (IMSExtract and GaussFit) and the .raw files containing your calibration spectra regardless of the polymer.

- Open a Windows command line (Shift + right click and select 'Open a command prompt here').

- $python IMS_Extract_2.0_PEG_PLA_POLYALA.py

- The script will detect which polymer it analyzes based on what's written in the header file.

- For each spectrum, the extracted arrival times are listed in the corresponding polymer folder and in the Output.log file. 
Only the reasonable fittings ($R^2 > 0.95$) are printed, the other ones are discarded and will have to be extracted by hand.


4. Additional notes

Currently supported species: [PEG + n Na]n+, [PLA + n Na]n+, [POLYALA + n H]n+ with n = 1+, 2+, 3+

Currently supported Mn: 600, 1000, 2000, 3350 g.mol-1 for PEG. The Mn of the two other polymers is less important.

The script detects the polymer type and size based on what's written in the header file. Conditions are applied depending on the charge state and approximate Mn. 
For instance, the program will not look for 3+ ions in PEG 600 as this polymer is too small to be triply charged.
!!! Again, specifying the polymer type and approximate Mn is mandatory !!!

The same ion can be found in two different spectrum (For instance, [PEG + 1 Na]1+ with 30 monomer units can be found in PEG 600 or PEG 1000 spectra). 
While the two ions' arrival times should be the same, the program might extract two slightly different arrival times. In this case, !!! check your data manually !!!

The first time you run IMSExtract, please check that the extracted ATDs match the arrival times from MassLynx.
We do not use the ATDs as provided from TWIMExtract because they might lead to a slight shift (max 1 ms) for Synapt G2-S and G2-Si instruments. 
The script instead reconstructs an arrival time scale based on the pusher frequency and offset which are read from the _extern.inf file.

Have fun calibrating!







