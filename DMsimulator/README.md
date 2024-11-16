This folder will eventually include the functions required to simulate a segmented DM.
It is currently organized as follows:

- The .yaml files contain the configuration parameters corresponding to the [TN] used for the simulation
All data generated during a simulation is saved in a folder with foldername corresponding to the configuration [TN]
If the same configuration file is used, data is READ (not COMPUTED) directly from the .fits files in the corresponding folder.
Thus, if you wish to change a parameter in a configuration file after running the code, you should either create a new .yaml file or delete the
corresponding folder before running the simulation

- geometry.py defines the hexagonal segments' geometry (e.g. masks, initial actuator location, segment disposition)

- main_script.py is where all the fun is happening
  
- read_configuration is a function to correctly read the parameters in the .yaml configuration file

- read_and_write_fits is a module containing functions to easily read/write masked arrays, block sparse matrices and lists in .fits format

- zernike_polynomials is a function to project the (you guessed it) zernike polynomials on a given mask
