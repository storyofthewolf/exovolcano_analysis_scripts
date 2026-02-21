"""
time_series.py - Orchestrator for exovolcano time series analysis.

Usage:
    python time_series.py ben2_vei7.yaml
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import config
import compute

print("\n!============= Running exovolcano time_series diagnostics =============!")

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

file_list = config.get_file_list()
if not file_list:
    raise SystemExit("No files found. Check config.")

exp_name    = config.get_experiment_name()
figures_dir = os.path.join(config.FIGURES_DIR, exp_name)
data_dir    = os.path.join(config.DATA_DIR,    exp_name)

os.makedirs(figures_dir, exist_ok=True)
os.makedirs(os.path.join(data_dir, 'scalar'),   exist_ok=True)
os.makedirs(os.path.join(data_dir, 'profiles'), exist_ok=True)

print(f"\nExperiment : {exp_name}")
print(f"Figures    : {figures_dir}")
print(f"Data       : {data_dir}\n")


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def save_scalar_csv(days, data_array, name, units=''):
    path = os.path.join(data_dir, 'scalar', f"{name}.csv")
    header = f"days_since_start,{name} [{units}]"
    np.savetxt(path, np.column_stack([days, data_array.values]),
               delimiter=',', header=header, comments='')
    print(f"  Saved scalar : {path}")


def save_profile_csv(days, data_array, name, pressure_1d, altitude_1d):
    """
    Row 0 = coordinate reference (P [Pa] and Z [m] per level).
    Subsequent rows = [day, val_lev0, val_lev1, ...].
    """
    path   = os.path.join(data_dir, 'profiles', f"{name}.csv")
    nlev   = len(pressure_1d)
    p_vals = np.asarray(pressure_1d)
    z_vals = np.asarray(altitude_1d) if altitude_1d is not None else np.full(nlev, np.nan)
    prof   = data_array.values   # (time, lev)

    lev_labels = [f"lev{k}" for k in range(nlev)]
    header = ('days_since_start,'
              + ','.join(f"P_{l}_Pa"   for l in lev_labels) + ','
              + ','.join(f"Z_{l}_m"    for l in lev_labels) + ','
              + ','.join(f"{name}_{l}" for l in lev_labels))

    coord_row = np.concatenate([[np.nan], p_vals, z_vals, np.full(nlev, np.nan)])
    rows = [np.concatenate([[day], p_vals, z_vals, prof[i, :]])
            for i, day in enumerate(days)]

    np.savetxt(path, np.vstack([coord_row, rows]),
               delimiter=',', header=header, comments='')
    print(f"  Saved profile: {path}")


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def plot_scalar_series(days, series_dict, title, ylabel, filename):
    fig, ax = plt.subplots(figsize=(12, 4))
    for label, da in series_dict.items():
        ax.plot(days, da.values, label=label)
    ax.set_xlabel('Days since start')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    path = os.path.join(figures_dir, filename)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved figure : {path}")


def plot_profile_hovmoller(days, data_array, pressure_1d, name,
                           filename, log_scale=False, units='',
):
    """
    Hovmoller diagram with log-pressure on the y-axis (surface at bottom,
    model top at top).

    log_scale     : bool  use LogNorm colormap (for trace species)
    Log-scale vars: vmax = data peak, vmin = vmax / 10^LOG_SCALE_DECADES
    Linear vars   : vmin/vmax clipped to 2nd-98th percentile
    """
    profile = data_array.values          # (time, lev)
    p_mb    = np.asarray(pressure_1d) / 100.0   # Pa -> mb

    fig, ax = plt.subplots(figsize=(12, 6))

    if log_scale:
        data_plot = np.where(profile > 0, profile, np.nan)
        pos       = data_plot[np.isfinite(data_plot) & (data_plot > 0)]
        decades   = LOG_SCALE_DECADES.get(name, 4)
        if pos.size > 0:
            _vmax = pos.max()
            _vmin = _vmax / 10**decades
        else:
            _vmax = _vmin = None
        norm  = mcolors.LogNorm(vmin=_vmin, vmax=_vmax) if (_vmin and _vmax) else None
        mesh  = ax.pcolormesh(days, p_mb, data_plot.T,
                              norm=norm, shading='auto', cmap='viridis')
    else:
        # Linear scale: clip at 2nd and 98th percentile to ignore outliers
        _vmin = float(np.nanpercentile(profile, 2))
        _vmax = float(np.nanpercentile(profile, 98))
        mesh = ax.pcolormesh(days, p_mb, profile.T,
                             vmin=_vmin, vmax=_vmax,
                             shading='auto', cmap='viridis')

    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label(f"{name} [{units}]" if units else name)

    # Log-pressure y-axis: surface (high P) at bottom, top (low P) at top
    ax.set_yscale('log')
    ax.set_ylim(p_mb.max(), p_mb.min())
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.4g'))
    ax.set_xlabel('Days since start')
    ax.set_ylabel('Pressure (mb)')
    ax.set_title(f'Global Mean {name} vs Time')

    plt.tight_layout()
    path = os.path.join(figures_dir, filename)
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Saved figure : {path}")


# Log-scale colormap settings for profile Hovmoller plots.
# Keys: variable name. Values: orders of magnitude below peak to show.
# Variables not listed here use linear scale with 2nd-98th percentile clipping.
LOG_SCALE_DECADES = {
    'SO2':     4,
    'H2SO4':   4,
    'Q':       4,
    'VOLCHZMD':4,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

ds, gw_values = compute.load_dataset(file_list)

with ds:
    days = compute.days_since_start(ds)

    print("Computing grid geometry...")
    geom = compute.compute_geometry(ds, gw_values)

    # Time-mean, area-mean pressure and altitude profiles for CSV output
    pressure_1d = geom['mid_p'].mean(dim=['time', 'lat', 'lon']).compute().values
    altitude_1d = geom['z_mid'].mean(dim=['time', 'lat', 'lon']).compute().values

    # -----------------------------------------------------------------------
    # Scalar time series
    # -----------------------------------------------------------------------
    print("\n--- Scalar time series ---")
    scalars = {}

    for var_cfg in config.SCALAR_VARS:
        name   = var_cfg['name']
        method = var_cfg['method']
        print(f"  Computing {name} ({method})...")
        result = compute.compute_scalar(ds, geom, name, method)
        if result is not None:
            scalars[name] = result
            units = ds[name].attrs.get('units', '') if name in ds else ''
            save_scalar_csv(days, result, name, units)

    # -----------------------------------------------------------------------
    # Profile time series
    # -----------------------------------------------------------------------
    print("\n--- Profile time series ---")
    profiles = {}

    for var_cfg in config.PROFILE_VARS:
        name = var_cfg['name']
        print(f"  Computing profile {name}...")
        result = compute.compute_profile(ds, geom, name)
        if result is not None:
            profiles[name] = result
            save_profile_csv(days, result, name, pressure_1d, altitude_1d)

    # Always include VOLCHZMD profile if variable exists
#    if 'VOLCHZMD' in ds.data_vars and 'VOLCHZMD' not in profiles:
#        print("  Computing profile VOLCHZMD...")
#        result = compute.compute_profile(ds, geom, 'VOLCHZMD')
#        if result is not None:
#            profiles['VOLCHZMD'] = result
#            save_profile_csv(days, result, 'VOLCHZMD', pressure_1d, altitude_1d)

    # -----------------------------------------------------------------------
    # Quicklook - scalar plots
    # -----------------------------------------------------------------------
    print("\n--- Quicklook plots ---")

    sulfur_vars = ['SO2', 'H2SO4', 'VOLCHZMD']
    sulfur = {k: scalars[k] for k in sulfur_vars if k in scalars}
    if sulfur:
        total = sum(sulfur.values())
        plot_scalar_series(days, dict(sulfur, Total=total),
                           title='Global Total Sulfur Budget',
                           ylabel='Mass (kg)',
                           filename='quicklook_sulfur_budget.png')

    for name in ['TS', 'TGCLDLWP', 'TMQ']:
        if name in scalars:
            units = ds[name].attrs.get('units', '') if name in ds else ''
            plot_scalar_series(days, {name: scalars[name]},
                               title=f'Global Mean {name}',
                               ylabel=units,
                               filename=f'quicklook_{name}.png')

    if 'Q' in scalars:
        plot_scalar_series(days, {'Q': scalars['Q']},
                           title='Global Total Water Vapor Mass',
                           ylabel='Mass (kg)',
                           filename='quicklook_Q_mass.png')

    # -----------------------------------------------------------------------
    # Quicklook - Hovmoller profile plots (log-pressure y-axis)
    # -----------------------------------------------------------------------
    for var_cfg in config.PROFILE_VARS + [{'name': n} for n in profiles if n not in [v['name'] for v in config.PROFILE_VARS]]:
        name = var_cfg['name']
        if name not in profiles:
            continue
        da    = profiles[name]
        units = ds[name].attrs.get('units', '') if name in ds else ''
        plot_profile_hovmoller(
            days, da, pressure_1d, name,
            filename=f'quicklook_profile_{name}.png',
            log_scale=(name in LOG_SCALE_DECADES),
            units=units,
        )

print(f"\nDone. Figures in '{figures_dir}', data in '{data_dir}'.")
print("\n!============= Exiting exovolcano time_series diagnostics =============!")
