import os
import glob

# --- !! USER CONFIGURATION !! ---

# 1. Path configuration
# Set your root directory here
ROOT_DIR = '/Users/wolfe/Desktop/projects/volcanos/runs'

# Set the filename pattern here (can be a single string or a list of strings)
FILE_PATTERN = [
    'eruption_test_land.cam.h1.0001-01-01-00000.nc',
    'eruption_test_land.cam.h1.0001-04-11-00000.nc',
    'eruption_test_land.cam.h1.0001-07-20-00000.nc',
    'eruption_test_land.cam.h1.0001-10-28-00000.nc',
    'eruption_test_land.cam.h1.0002-02-05-00000.nc',
    'eruption_test_land.cam.h1.0002-05-16-00000.nc',
    'eruption_test_land.cam.h1.0002-08-24-00000.nc',
    'eruption_test_land.cam.h1.0002-12-02-00000.nc',
    'eruption_test_land.cam.h1.0003-03-12-00000.nc',
    'eruption_test_land.cam.h1.0003-06-20-00000.nc',
    'eruption_test_land.cam.h1.0003-09-28-00000.nc'
]


# 2. Directory to save plot
OUTPUT_DIR = 'figures'

# --- Physical Constants ---
#G_CONST = 9.80665      # Gravitational acceleration (m/s^2)
#R_AIR = 287.058        # Specific gas constant for dry air (J/kg*K)
#R_EARTH = 6.371e6      # Mean radius of Earth (m)

# t1e pure CO2
G_CONST = 9.80665*0.93      # Gravitational acceleration (m/s^2)
R_AIR = 188.965172522727 # pure CO2 atmosphere
R_EARTH = 6.371e6*0.91      # Mean radius of Earth (m)

# --- Logic to generate file list ---
def get_file_list():
    """Discovers and returns a sorted list of unique files based on configuration."""
    if isinstance(FILE_PATTERN, str):
        patterns = [FILE_PATTERN]
    else:
        patterns = FILE_PATTERN

    file_list = []
    for pattern in patterns:
        full_path = os.path.join(ROOT_DIR, pattern)
        expanded_path = os.path.expanduser(full_path)
        print(f"Searching for files with pattern: {expanded_path}")
        file_list.extend(glob.glob(expanded_path))

    # Remove duplicates (if patterns overlap) and sort
    file_list = sorted(list(set(file_list)))
    
    if not file_list:
        print("\n" + "!"*50)
        print("ERROR: No NetCDF files found!")
        print(f"Patterns searched in: {ROOT_DIR}")
        print(f"Current Working Directory: {os.getcwd()}")
        print("Please check if ROOT_DIR and FILE_PATTERN in config.py are correct.")
        print("!"*50 + "\n")
        return []
        
    return file_list
