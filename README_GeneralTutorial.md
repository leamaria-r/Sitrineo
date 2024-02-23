# Sitrineo Analysis M2 2024  

The goal of this project is to use the Sitrineo (Silicon Tracker with International Education Objective) tracker and measure the energy spectrum emitted by a Sr-90 radioactive source (beta emitter). Data are analysed using the TAF software, based on C++ and ROOT. This is a master students project made by Mattéo Maushart and Léa-Maria Rabour. For more details, see https://github.com/jeromebaudot/taf

## Repositories

### code
contains the code of TAF software, separated in different directories :
- **include** with all the .h files included in the C++ codes in source
- **src** with all the source code
- **macros** with some C macros

README explains the modifications made to the existing code.

### Config
contains the config used to analyse data with TAF

### Results
contains the results of the TAF analysis. Each run has its own directory

### Runs
contains the data acquired by Sitrineo. Each run has its own directory

### Sitrineo Physics
contains the code to simulate the activity of the source and the energy loss in the tracker

### Plots
contains the main plots obtained with the analysis

## How to run the analysis ?  

The requirements to execute the code are :  
    -root : version 6.10/04
    -gcc : version 7.5

The computer's operating system is Ubuntu version 7.5 (3ubuntu1~18.04)
 
### How to use TAF software?

Open a terminal and cd in the folder where taf is installed (e.g. cd home/dphe1/physics/taf).  
Run taf with the command : taf -run [RRRR] -cfg [config path] (e.g. taf -run 1110 -cfg ./config/sitrineo-m2-B.cfg).  
NB : to open the graphic interface, one need to add -guiw at the end of the command (e.g. taf -run 1110 -cfg ./config/sitrineo-m2-B.cfg -guiw).  
If running taf without graphical interface, one need to call the function with the command : gTAF->GetRaw()->[function(parameter)] (e.g. gTAF->GetRaw()->SitrineoAlign(100000)).  
The plots are displayed on screen and saved in taf/Results/ .  

### Name of the runs

We name the runs with up to 4 numbers (RRRR) using the following convention :   
1st number (R...) is equal to 1 if there is a magnetic field B, 0 without B   
2nd number (.R..) is equal to the diameter of the collimator in mm (0 for runs without collimator)   
3rd number (..R.) is equal to 1 if there is a source, 0 without the source   
4th number (...R) is the number of runs with such config before this run (i.e. 0 for the first run with a given config, 1 for the second etc.)   
Like for usual numbers, 0 on the left are not written so run 0010 would be run 10.   


For example, the first run with source, magnetic field and 1mm collimator will be named 1110. The second one would be 1111, and so on.   
NB: Event 0 cannot exist (otherwise taf isn't happy...) so a run without source, magnetic field nor collimator would exceptionally be named 1, the second one 2, etc.   

### Functions that are used

* SitrineoAlign() is used to align the planes of the tracker. We automatised the procedure so the config file is now updated with the new positions that are computed. It takes as an argument a number of events (e.g. 100000). This function must be used on a config where B=0.0 (no magnetic field) and a run without magnetic field, with a source and a collimator (e.g. run 115 or 116)
* SitrineoCumul() is used to fit the tracks and plot the momentum spectrum. It takes as an argument a number of events (e.g. 100000). This function must be used on a config file where B!=0 (here B=0.2) and after the planes were aligned. Easiest way to do that is to update the value of BFieldMagnitude (line 54 of the config file, value must be a float) in the config file that has been used by SitrineoAlign() before. 
The function also requires a run with magnetic field, a source and a collimator (e.g. run 1110 or 1111)
Warning: not displayed in graphical interface (for now)
* DisplayCumulatedRawData2D() (or "Create noisy pixel map" if using graphical interface) is used to create a map of the noisy pixels. It must be used with a run without magnetic field nor source (hence no collimator either), e.g. 1 or 2, and a high number of events (e.g. 200000). This map can then be used to mask hot pixels.


## Data acquisition with Sitrineo 

First, one need to turn on the power supply and turn on the Sitrineo DAQBoard. After waiting a few seconds, connect your computer to the DAQBoard by entering the following command in the terminal : ssh root@192.168.1.12   
Then enter the following commands to start data acquisition :  
sitrireset  
sitriconf   
sitristart  
daqSoC_v3 -RunNumber 1210 -NumEventsToRead 100000 -triggerSW yes -SteptriggerMonitor 2000 -JtagInit work -Delay 1 -DataSaveLocal yes  
(for run 1210, otherwise change run number) 

After that, you can open a new terminal and go to your own repository to copy the files there with :  
scp -r root@192.168.1.12:DataStore/run1213 1213   
(for run 1213, otherwise change run number)
