This folder will eventually include the functions required to simulate a segmented DM.
It is currently organized as follows:

- The .yaml files contain the configuration parameters corresponding to the [TN] used for the simulation
All data generated during a simulation is saved in a folder with foldername corresponding to the configuration [TN]
If the same configuration file is used, data is READ (not COMPUTED) directly from the .fits files in the corresponding folder
This means that, if you wish to change a parameter in a configuration file, you should either create a new .yaml file or delete the
corresponding folder before running the simulation

- geometry.py defines the hexagonal segments' geometry (e.g. masks, interaction matrix)

- main.py is where all the fun is happening

- all other files are quite self explanatory (I'm sure you can figure out what read_config.py does)