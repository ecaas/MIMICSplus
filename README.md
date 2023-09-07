# MIMICS+
This repository contains the source code to run the soil decomposition model MIMICS+. The [test](./test) folder contains inputdata for a test case (Hurdal) with which you should be able to run the model following the quick start guide below.

## Quick start  
1. Clone repository.  
2. To compile you need netcdf-fortran package. How to create a new conda environment, and download the package: 
 - "conda create --name myenv"
 - "conda activate myenv"
 - "conda install -c conda-forge netcdf-fortran"   
You can also use the environment.yml file. 
3. cd into the repository and open the Makefile. Edit the paths to /lib and /include so that they point to where the netcdf-fortran package is installed. Save and exit.  
4. Compile the code by typing "make" (make sure you have make installed) in the command line. Hopefully this will work without errors. 
5. If point 4. worked, you can open the "run.bash" file and set name of simulation, how long the spinup period should be, and if you want to run historical or just spinup. If you want to use other input files than those in the "test" directory, you need to change the paths here.  
6. Start simulation by typing "bash run.bash" in the command line. Output files should appear in "results/name_of_simulation/".
7. The simulation length (nsteps) and output frequency (write_hour) can be adjusted by the calls to the main subroutine in [main.F90](https://github.com/ecaas/MIMICSplus/blob/93e45deb9f7b3e9a4edfcabc985b0a01cc42f8da/src/main.f90#L138C1-L169C8) 
