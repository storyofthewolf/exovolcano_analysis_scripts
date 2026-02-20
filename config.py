"""
config.py - Experiment configuration loader for exovolcano analysis.

Reads experiment parameters from a YAML file. The config file path defaults
to 'experiment.yaml' in the current directory, or can be specified via the
CONFIG environment variable or as the first command-line argument.

Usage:
    python time_series.py                        # uses experiment.yaml
    python time_series.py ben2_vei7.yaml         # uses named file
    CONFIG=ben2_vei7.yaml python time_series.py  # via environment variable
"""

import os
import sys
import glob
import yaml


# ---------------------------------------------------------------------------
# Load YAML config file
# ---------------------------------------------------------------------------
EXPERIMENTS_DIR = "experiments"


def _find_config_file():
    """Resolve config file path from CLI arg, env var, or default."""
    if len(sys.argv) > 1 and sys.argv[1].endswith('.yaml'):
        name = sys.argv[1]
        # Accept full path, or bare name resolved into experiments/
        if os.path.exists(name):
            return name
        return os.path.join(EXPERIMENTS_DIR, name)
    if 'CONFIG' in os.environ:
        return os.environ['CONFIG']
    return os.path.join(EXPERIMENTS_DIR, 'experiment.yaml')


def _load_config(path):
    """Load and return the YAML config as a dict."""
    if not os.path.exists(path):
        print(f"ERROR: Config file not found: '{path}'")
        sys.exit(1)
    with open(path) as f:
        cfg = yaml.safe_load(f)
    # Normalize file_pattern to always be a list
    if isinstance(cfg.get('file_pattern'), str):
        cfg['file_pattern'] = [cfg['file_pattern']]
    return cfg


_cfg = _load_config(_find_config_file())


# ---------------------------------------------------------------------------
# Public constants (same names as before — time_series.py unchanged)
# ---------------------------------------------------------------------------

ROOT_DIR    = _cfg['root_dir']
FOLDER      = _cfg['folder']
FILE_PATTERN = _cfg['file_pattern']
G_CONST     = _cfg['g_const']
R_AIR       = _cfg['r_air']
R_EARTH     = _cfg['r_earth']
OUTPUT_DIR  = _cfg.get('output_dir', 'figures')


# ---------------------------------------------------------------------------
# Public helper functions (identical API to original config.py)
# ---------------------------------------------------------------------------

def get_file_list():
    """Discovers and returns a sorted list of files based on configuration."""
    file_list = []
    for pattern in FILE_PATTERN:
        full_path = os.path.join(ROOT_DIR, FOLDER, pattern)
        expanded = os.path.expanduser(full_path)
        print(f"Searching: {expanded}")
        file_list.extend(glob.glob(expanded))

    file_list = sorted(set(file_list))

    if not file_list:
        print("\n" + "!" * 50)
        print("ERROR: No NetCDF files found!")
        print(f"  root_dir:     {ROOT_DIR}")
        print(f"  folder:       {FOLDER}")
        print(f"  file_pattern: {FILE_PATTERN}")
        print("!" * 50 + "\n")

    return file_list


def get_experiment_name():
    """Returns the common filename prefix shared by all files in FILE_PATTERN."""
    prefixes = {f.split('.')[0] for f in FILE_PATTERN}

    if len(prefixes) > 1:
        print(f"ERROR: Multiple file prefixes detected: {prefixes}")
        sys.exit(1)
    if len(prefixes) == 0:
        print("ERROR: file_pattern is empty.")
        sys.exit(1)

    name = list(prefixes)[0]
    print(f"Experiment name: {name}")
    return name
