# ParaView_NetCDF_Particles_Reader_Plugin
A simple, user-friendly Python plugin for ParaView that reads particles in NetCDF format. The plugin was developed and tested with ParaView 5.12 and Python 3.10.

Please read the tutorial before using the plugin (NetCDF_ParticlesReader_TUTORIAL.pdf).

We provided the example files shown in the tutorial for self-study:  
- TERIFIC_drifters5_Github.nc - extracted and processed from a data set courtesy of Prof. Eleanor Frajka-Williams (Institut für Meereskunde, Universität Hamburg)
- floats_GULF_2007_1510STRIDE.nc - extracted and processed from a data set courtesy of Dr. Alexa Griesel (Institut für Meereskunde, Universität Hamburg)
- radar_7000TS_ALTITUDE.nc - extracted and processed from an EUREC4A data set https://observations.ipsl.fr/thredds/dodsC/EUREC4A/PRODUCTS/TRACKS/EUREC4A_tracks_HALO_v1.0.nc,
                             paper: Konow et al. (2021). EUREC4A’s HALO. Earth System Science Data, 13(12), 5545–5563. https://doi.org/10.5194/essd-13-5545-2021

My work Python - ParaView configuration on Windows 10: 
- First I installed the neCDF-4 Python module in a Python 3.10.13 environment in Anaconda
- The ParaView version I used is 5.12, which also has Python 3.10
- When working, open an Anaconda prompt
-   then activate the relevant environment, in my case :
    conda activate PYTHON310
-   and set up the PYTHONHOME and PYTHONPATH variables, for ParaView to work with this environment and recognize the necessary libraries:
    set PYTHONHOME=my_path_to\anaconda3\envs\PYTHON310
    
    set PYTHONPATH=my_path_to\anaconda3\envs\PYTHON310\lib
    
-   and finally, from this same Anaconda prompt, start ParaView:    
    my_path_to\ParaView-5.12.0\bin\paraview.exe

