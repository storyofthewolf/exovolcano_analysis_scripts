"""
zonal_plots.py - Zonal mean contour plot for exovolcano analysis.

Public API:
    plot_zonal_mean(lat, pressure_1d, data_2d, name, units,
                    actual_day, figures_dir, filename, log_scale=False)

LOG_SCALE_DECADES is also exported; run_time_series.py imports it from here
so the definition lives in one place.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# Variables that use LogNorm colormaps and their dynamic range in decades.
# All other variables use linear scale with 2nd-98th percentile clipping.
LOG_SCALE_DECADES = {
    'SO2':      4,
    'H2SO4':    4,
    'Q':        4,
    'VOLCHZMD': 4,
}


def plot_zonal_mean(lat, pressure_1d, data_2d, name, units,
                    actual_day, figures_dir, filename, log_scale=False):
    """
    Zonal mean contour plot. x=latitude, y=log-pressure (surface at bottom,
    top at top), color=variable value.

    Parameters
    ----------
    lat          : np.ndarray (nlat,)   latitude [degrees]
    pressure_1d  : np.ndarray (nlev,)   pressure [Pa], time+area mean
    data_2d      : np.ndarray (lev, lat)
    name         : str   variable name (used in title and colorbar)
    units        : str   variable units (used in colorbar label)
    actual_day   : float day of the selected timestep
    figures_dir  : str   output directory
    filename     : str   output PNG filename
    log_scale    : bool  use LogNorm anchored at data peak
    """
    p_mb = np.asarray(pressure_1d) / 100.0   # Pa -> mb

    fig, ax = plt.subplots(figsize=(12, 6))

    if log_scale:
        plot_data = np.where(data_2d > 0, data_2d, np.nan)
        pos       = plot_data[np.isfinite(plot_data) & (plot_data > 0)]
        decades   = LOG_SCALE_DECADES.get(name, 4)
        if pos.size > 0:
            _vmax = pos.max()
            _vmin = _vmax / 10**decades
        else:
            _vmax = _vmin = None
        norm = mcolors.LogNorm(vmin=_vmin, vmax=_vmax) if (_vmin and _vmax) else None
        mesh = ax.pcolormesh(lat, p_mb, plot_data,
                             norm=norm, shading='auto', cmap='viridis')
    else:
        _vmin = float(np.nanpercentile(data_2d, 2))
        _vmax = float(np.nanpercentile(data_2d, 98))
        mesh = ax.pcolormesh(lat, p_mb, data_2d,
                             vmin=_vmin, vmax=_vmax,
                             shading='auto', cmap='viridis')

    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label(f"{name} [{units}]" if units else name)

    ax.set_yscale('log')
    ax.set_ylim(p_mb.max(), p_mb.min())
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.4g'))
    ax.set_xlim(lat.min(), lat.max())
    ax.set_xlabel('Latitude (degrees)')
    ax.set_ylabel('Pressure (mb)')
    ax.set_title(f'Zonal Mean {name} — Day {actual_day:.2f}')

    plt.tight_layout()
    path = os.path.join(figures_dir, filename)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved figure : {path}")
