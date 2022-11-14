# soil_decomp
## Quick start  
1. Clone repository.  
2. To compile you need netcdf-fortran package. How to create a new conda environment, and download the package: 
 - "conda create --name myenv"
 - "conda activate myenv"
 - "conda install -c conda-forge netcdf-fortran"   
You can also use the environment.yml file (NB: Not fully tested that the file works!). 
3. cd into the repository and open the Makefile. Edit the paths to /lib and /include so that they point to where the netcdf-fortran package is installed. Save and exit.  
4. Compile the code by typing "make" in the command line. Hopefully this will work without errors. 
5. If point 4. worked, you can open the "run.bash" file and set name of simulation, how long the spinup period should be, and if you want to run historical or just spinup. If you want to use other input files than those in the "test" directory, you need to change the paths here.  
6. Start simulation by typing "bash run.bash" in the command line. Output files should appear in "results/name_of_simulation/".
