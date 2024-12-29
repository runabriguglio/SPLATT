"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024
    - Matteo Menessini : adapted for SPLATT in 12/2024

Description
-----------
Builds the paths folder tree to be used for splatt

"""
import os
import sys

try:
    BASE_PATH = os.getenv('SPLATTBASEPATH')
except KeyError as exc:
    raise KeyError("Environment variable not found! Define the SPLATTBASEPATH env variable that points to '~/SPLATT/splattsw") from exc

# BASE_DATA_PAT      = os.path.join(BASE_PATH, 'data')
# WEBDAQ_PATH        = os.path.join(BASE_PATH, 'DMSimulator')
DM_SIMULATOR_PATH    = os.path.join(BASE_PATH, 'DMsimulator')
SPLATT_SOFTWARE_PATH = os.path.join(BASE_PATH, 'splattsw')
SPLATT_DEVICES_PATH  = os.path.join(SPLATT_SOFTWARE_PATH, 'devices')

paths = [DM_SIMULATOR_PATH, SPLATT_SOFTWARE_PATH, SPLATT_DEVICES_PATH]
for p in paths:
    if not os.path.exists(p):
        os.mkdir(p)
    
    # Add to path
    sys.path.append(p)