"""
aod_plots.py - AOD plotting functions for exovolcano analysis.

All functions accept a figures_dir argument and write a PNG there.
Callers are expected to pass the aod/ subdirectory
(figures/<exp_name>/aod/) so AOD figures are kept separate from the
top-level quicklook plots.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_aod_timeseries(days, aod_global_mean, wavelength_label, figures_dir, filename):
    """
    Global-mean AOD time series.

    Parameters
    ----------
    days             : np.ndarray (time,)  days since start
    aod_global_mean  : np.ndarray (time,)  area-weighted global mean AOD
    wavelength_label : str                 e.g. '550 nm (band)' or '0.55 µm (Mie)'
    figures_dir      : str                 output directory
    filename         : str                 output filename (PNG)
    """
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(days, aod_global_mean)
    ax.set_xlabel('Days since start')
    ax.set_ylabel('AOD')
    ax.set_title(f'Global Mean AOD at {wavelength_label}')
    ax.grid(True)
    plt.tight_layout()
    path = os.path.join(figures_dir, filename)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved figure : {path}")


def plot_aod_zonal_hovmoller(days, aod_zonal_mean, lat, figures_dir, filename):
    """
    Zonal-mean AOD Hovmoller diagram (time x latitude).

    Parameters
    ----------
    days           : np.ndarray (time,)       days since start
    aod_zonal_mean : np.ndarray (time, lat)   zonal-mean AOD
    lat            : np.ndarray (nlat,)        latitude in degrees
    figures_dir    : str                       output directory
    filename       : str                       output filename (PNG)

    Notes
    -----
    # To use equal-area spacing replace lat with np.sin(np.deg2rad(lat))
    # and update ylabel to 'sin(latitude)'
    """
    _vmin = float(np.nanpercentile(aod_zonal_mean, 2))
    _vmax = float(np.nanpercentile(aod_zonal_mean, 98))

    fig, ax = plt.subplots(figsize=(12, 6))
    mesh = ax.pcolormesh(days, lat, aod_zonal_mean.T,
                         vmin=_vmin, vmax=_vmax,
                         shading='auto', cmap='viridis')
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label('AOD')
    ax.set_xlabel('Days since start')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title('Zonal Mean AOD vs Time')
    plt.tight_layout()
    path = os.path.join(figures_dir, filename)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved figure : {path}")
