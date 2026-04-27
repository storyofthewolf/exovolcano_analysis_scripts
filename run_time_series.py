#!/usr/bin/env python

"""
run_time_series.py - Orchestrator for exovolcano time series analysis.

Usage:
    python run_time_series.py ben2_vei7.yaml
"""

import sys
import time as _time
import datetime as _datetime

_t_start_wall = _time.perf_counter()
_t_start_dt   = _datetime.datetime.now()

# Parse all custom flags before config.py sees sys.argv
_timing_enabled = '--time' in sys.argv
if _timing_enabled:
    sys.argv.remove('--time')

_SKIP_FLAGS = {'--no-zonal', '--no-aod', '--no-profiles', '--no-scalars', '--no-plots'}
_skips = _SKIP_FLAGS & set(sys.argv)
for _f in _skips:
    sys.argv.remove(_f)

_n_workers = 8
if '--nthreads' in sys.argv:
    _i = sys.argv.index('--nthreads')
    _n_workers = int(sys.argv[_i + 1])
    sys.argv[_i:_i + 2] = []

_skip_scalars  = '--no-scalars'  in _skips
_skip_profiles = '--no-profiles' in _skips
_skip_plots    = '--no-plots'    in _skips
_skip_aod      = '--no-aod'      in _skips
_skip_zonal    = '--no-zonal'    in _skips

_timings: dict = {}   # label -> elapsed seconds (insertion order = print order)

def _tick(label):
    """Start a named timer. Returns a stop callable."""
    if not _timing_enabled:
        return lambda: None
    t0 = _time.perf_counter()
    def _tock():
        _timings[label] = _timings.get(label, 0.0) + (_time.perf_counter() - t0)
    return _tock

if len(sys.argv) > 1 and sys.argv[1] in ('-h', '--help'):
    print("""
Usage:
    python run_time_series.py <experiment>.yaml
    python run_time_series.py experiments/<experiment>.yaml
    CONFIG=experiments/<experiment>.yaml python run_time_series.py

  The YAML file is looked up in the experiments/ directory unless a path
  separator is present.  A fully-annotated template is at:
    experiments/template.yaml

Required YAML keys
------------------
  root_dir      Absolute path to the top-level CESM archive directory.
  folder        Subdirectory under root_dir containing the NetCDF files
                (e.g. 'run_name/atm/hist').
  file_pattern  List of CAM h1 history filenames (basenames).  All files
                must share the same prefix; that prefix becomes the
                experiment name used for output subdirectories.
  g_const       Surface gravity [m/s²].
  r_air         Specific gas constant of the atmosphere [J/kg/K].
  r_earth       Mean planet radius [m].

Optional YAML keys (defaults in parentheses)
--------------------------------------------
  figures_dir   Root directory for PNG output  ('figures').
  data_dir      Root directory for CSV output  ('data').
  optics_file   Path to volc_pw1975_n68_r1.0um_mie.nc.  Set this to enable AOD
                calculations; omit or set to null to skip AOD entirely.
  volc_reff     Effective particle radius [µm] for Kext lookup  (1.0).
  rho_aerosol   Bulk aerosol density [g/cm³], used by Mie path  (1.84).
  mie_wavelength_um          Wavelength [µm] for optional single-wavelength
                             Mie AOD calculation.  Requires miepython.
  mie_refractive_index_real  Real part of complex refractive index  (1.43).
  mie_refractive_index_imag  Imaginary part, positive convention  (0.0).

scalar_vars / profile_vars sections
------------------------------------
  Each entry is a list item with 'name' and 'method' keys:

    scalar_vars:
      - name: SO2
        method: mass_integral     # global total mass [kg]
      - name: TS
        method: area_mean         # Gaussian-weighted global mean
      - name: VOLCHZMD
        method: volume_integral   # density-weighted total mass [kg]

    profile_vars:
      - name: T
        method: area_mean         # global mean profile (time × lev)

  method options:  mass_integral | volume_integral | area_mean

Runtime flags (can be combined)
--------------------------------
  --time           Print a per-section timing summary at the end.
  --nthreads N     Use N dask threads (default: 8). Increase on an allocated
                   compute node; keep low (4-8) on a shared login node.
  --no-scalars     Skip scalar time series (computation + CSV output).
  --no-profiles    Skip profile time series (computation + CSV output).
  --no-plots       Skip all figure output (CSV data is still written).
  --no-aod         Skip AOD calculation.
  --no-zonal       Skip zonal mean snapshots (expensive on large datasets).

Example
-------
    python run_time_series.py t1d_vei7.yaml
    python run_time_series.py t1d_vei7.yaml --nthreads 32 --no-zonal --time
""")
    sys.exit(0)

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import xarray as xr
import dask
import dask.config
import pandas as pd

dask.config.set(scheduler='threads', num_workers=_n_workers)

import config
import compute
import optics
import aod_plots
import zonal_plots
from zonal_plots import LOG_SCALE_DECADES

print("\n!============= Running exovolcano time_series diagnostics =============!")
print(f"Dask scheduler  : threaded  ({_n_workers} threads)")
if _skips:
    print(f"Skipping        : {' '.join(f.replace('--no-', '') for f in sorted(_skips))}")

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

_stop = _tick("Setup")
file_list = config.get_file_list()
if not file_list:
    raise SystemExit("No files found. Check config.")

exp_name    = config.get_experiment_name()
figures_dir = os.path.join(config.FIGURES_DIR, exp_name)
data_dir    = os.path.join(config.DATA_DIR,    exp_name)

os.makedirs(figures_dir, exist_ok=True)
os.makedirs(os.path.join(data_dir,    'scalar'),   exist_ok=True)
os.makedirs(os.path.join(data_dir,    'profiles'), exist_ok=True)
os.makedirs(os.path.join(data_dir,    'aod'),      exist_ok=True)
os.makedirs(os.path.join(figures_dir, 'aod'),      exist_ok=True)
_stop()

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
    Vertical coordinate info (P and Z per level) is written as comment lines
    before the column header so it is preserved without NaN placeholders.
    Data rows: days_since_start, val_lev0, val_lev1, ...
    """
    path   = os.path.join(data_dir, 'profiles', f"{name}.csv")
    nlev   = len(pressure_1d)
    p_vals = np.asarray(pressure_1d)
    z_vals = np.asarray(altitude_1d) if altitude_1d is not None else np.full(nlev, np.nan)
    prof   = data_array.values   # (time, lev)

    lev_labels = [f"lev{k}" for k in range(nlev)]
    header = ('days_since_start,'
              + ','.join(f"{name}_{l}" for l in lev_labels))

    with open(path, 'w') as f:
        f.write('# pressure_Pa: ' + ','.join(f"{v:.6e}" for v in p_vals) + '\n')
        f.write('# altitude_m: '  + ','.join(f"{v:.6e}" for v in z_vals) + '\n')
        f.write(header + '\n')
        np.savetxt(f, np.column_stack([days, prof]), delimiter=',')
    print(f"  Saved profile: {path}")


def save_aod_scalar_csv(days, aod_global, tag):
    """Global-mean AOD time series.  Two columns: days_since_start, AOD."""
    path   = os.path.join(data_dir, 'aod', f"aod_{tag}.csv")
    header = f"days_since_start,AOD_{tag}"
    np.savetxt(path, np.column_stack([days, aod_global]),
               delimiter=',', header=header, comments='')
    print(f"  Saved AOD scalar : {path}")


def save_aod_zonal_csv(days, aod_zonal, lat, tag):
    """Zonal-mean AOD.  Header row lists lat values; each data row is one timestep."""
    path      = os.path.join(data_dir, 'aod', f"aod_zonal_{tag}.csv")
    lat_labels = [f"{v:.4f}" for v in lat]
    header    = 'days_since_start,' + ','.join(lat_labels)
    np.savetxt(path, np.column_stack([days, aod_zonal]),
               delimiter=',', header=header, comments='')
    print(f"  Saved AOD zonal  : {path}")


def save_zonal_csv(data_2d, pressure_1d, lat_vals, name, actual_day):
    """
    Zonal mean snapshot. Rows=pressure levels, columns=latitudes.
    First row: lat header. First column: pressure in mb.
    """
    path = os.path.join(data_dir, 'zonal', f"{name}_day{actual_day:07.2f}.csv")
    p_mb = pressure_1d / 100.0
    cols = {f"{v:.4f}": data_2d[:, i] for i, v in enumerate(lat_vals)}
    df   = pd.DataFrame({'pressure_mb': p_mb, **cols})
    df.to_csv(path, index=False, float_format='%.6e')
    print(f"  Saved zonal  : {path}")


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



# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_stop = _tick("Load dataset")
ds, gw_values = compute.load_dataset(file_list)
_stop()

with ds:
    days = compute.days_since_start(ds)

    print("Computing grid geometry...")
    _stop = _tick("Geometry")
    geom = compute.compute_geometry(ds, gw_values)

    # Area-mean pressure and altitude profiles for CSV/plot labeling.
    # Use first timestep only — time variation is negligible for axis labels
    # and this avoids a full-dataset scan over all timesteps.
    _p1d, _z1d = dask.compute(
        geom['mid_p'].isel(time=0).mean(dim=['lat', 'lon']),
        geom['z_mid'].isel(time=0).mean(dim=['lat', 'lon']),
    )
    pressure_1d = _p1d.values
    altitude_1d = _z1d.values
    _stop()

    # -----------------------------------------------------------------------
    # Scalar time series
    # -----------------------------------------------------------------------
    print("\n--- Scalar time series ---")
    scalars = {}

    if _skip_scalars:
        print("  Skipped (--no-scalars).")
    else:
        _stop = _tick("Scalar time series")
        _scalar_lazy  = {}
        _scalar_units = {}
        for var_cfg in config.SCALAR_VARS:
            name   = var_cfg['name']
            method = var_cfg['method']
            print(f"  Computing {name} ({method})...")
            result = compute.compute_scalar(ds, geom, name, method)
            if result is not None:
                _scalar_lazy[name]  = result
                _scalar_units[name] = ds[name].attrs.get('units', '') if name in ds else ''

        if _scalar_lazy:
            _computed = dask.compute(*_scalar_lazy.values())
            scalars = dict(zip(_scalar_lazy.keys(), _computed))
            for name, result in scalars.items():
                save_scalar_csv(days, result, name, _scalar_units[name])
        _stop()

    # -----------------------------------------------------------------------
    # Profile time series
    # -----------------------------------------------------------------------
    print("\n--- Profile time series ---")
    profiles = {}

    if _skip_profiles:
        print("  Skipped (--no-profiles).")
    else:
        _stop = _tick("Profile time series")
        _profile_lazy = {}
        for var_cfg in config.PROFILE_VARS:
            name = var_cfg['name']
            print(f"  Computing profile {name}...")
            result = compute.compute_profile(ds, geom, name)
            if result is not None:
                _profile_lazy[name] = result

        if _profile_lazy:
            _computed = dask.compute(*_profile_lazy.values())
            profiles = dict(zip(_profile_lazy.keys(), _computed))
            for name, result in profiles.items():
                save_profile_csv(days, result, name, pressure_1d, altitude_1d)
        _stop()

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

    if _skip_plots:
        print("  Skipped (--no-plots).")
    else:
        _stop = _tick("Quicklook plots")
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

        # Hovmoller profile plots (log-pressure y-axis)
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
        _stop()

    # -----------------------------------------------------------------------
    # AOD
    # -----------------------------------------------------------------------
    print("\n--- AOD ---")

    if _skip_aod:
        print("  Skipped (--no-aod).")
    elif config.OPTICS_FILE is None:
        print("  WARNING: 'optics_file' not set in config. Skipping AOD calculations.")
    elif 'VOLCHZMD' not in ds.data_vars:
        print("  WARNING: 'VOLCHZMD' not found in dataset. Skipping AOD calculations.")
    else:
        _stop = _tick("AOD")
        band_optics    = optics.load_band_optics(config.OPTICS_FILE)
        volchzmd_vals  = ds['VOLCHZMD'].values   # (time, lev, lat, lon) [g/cm³]
        dz_vals        = geom['dz'].values        # (time, lev, lat, lon) [m]
        lat            = ds['lat'].values

        aod_figures_dir = os.path.join(figures_dir, 'aod')

        def _aod_diagnostics(kext, label, tag):
            """Compute global-mean and zonal-mean AOD, save CSVs and plots."""
            aod_2d = optics.compute_aod(volchzmd_vals, dz_vals, kext)  # (time, lat, lon)

            # Area-weighted global mean: weight by Gaussian weights over lat,
            # then simple mean over lon (all longitudes are equal-weight).
            aod_da  = xr.DataArray(aod_2d,
                                   dims=['time', 'lat', 'lon'],
                                   coords={'lat': ds['lat'], 'lon': ds['lon']})
            weights = xr.DataArray(gw_values, coords={'lat': ds['lat']}, dims=['lat'])
            aod_global = aod_da.weighted(weights).mean(dim=['lat', 'lon']).values

            # Zonal mean: simple mean over lon
            aod_zonal = aod_2d.mean(axis=2)   # (time, lat)

            print(f"  {label}: peak global-mean AOD = {aod_global.max():.4f}")

            save_aod_scalar_csv(days, aod_global, tag)
            save_aod_zonal_csv(days, aod_zonal, lat, tag)

            aod_plots.plot_aod_timeseries(
                days, aod_global, label, aod_figures_dir,
                filename=f'aod_{tag}_timeseries.png',
            )
            aod_plots.plot_aod_zonal_hovmoller(
                days, aod_zonal, lat, aod_figures_dir,
                filename=f'aod_{tag}_zonal_hovmoller.png',
            )

        # Band-interpolated Kext at 550 nm
        i_wave    = optics.select_band_550nm(band_optics['wvn_centers'])
        kext_550  = optics.interpolate_kext(band_optics, i_wave, config.VOLC_REFF)
        print(f"  Band 550 nm index={i_wave}, "
              f"center={band_optics['wvn_centers'][i_wave]:.1f} cm⁻¹, "
              f"Kext={kext_550:.4f} cm²/g  (reff={config.VOLC_REFF} µm)")
        _aod_diagnostics(kext_550, '550 nm (band)', '550nm_band')

        # Optional single-wavelength Mie Kext
        if config.MIE_WAVE_UM is not None:
            mie_wave = float(config.MIE_WAVE_UM)
            kext_mie = optics.mie_kext(
                mie_wave, config.VOLC_REFF,
                config.MIE_N_REAL, config.MIE_N_IMAG,
                config.RHO_AEROSOL,
            )
            label_mie = f'{mie_wave:.3f} µm (Mie)'
            print(f"  Mie {label_mie}: Kext={kext_mie:.4f} cm²/g")
            tag_mie = f'{mie_wave:.3f}um_mie'.replace('.', 'p')
            _aod_diagnostics(kext_mie, label_mie, tag_mie)
        _stop()

    # -----------------------------------------------------------------------
    # Zonal mean snapshots
    # -----------------------------------------------------------------------
    print("\n--- Zonal mean snapshots ---")

    if _skip_zonal:
        print("  Skipped (--no-zonal).")
    elif config.ZONAL_VARS and config.ZONAL_PERIODS:
        os.makedirs(os.path.join(data_dir,    'zonal'), exist_ok=True)
        os.makedirs(os.path.join(figures_dir, 'zonal'), exist_ok=True)

        lat_vals       = ds['lat'].values
        zonal_fig_dir  = os.path.join(figures_dir, 'zonal')

        _stop = _tick("Zonal mean snapshots")

        # Collect lazy lon-means, then single dask pass for all zonal vars.
        _zonal_valid = [v for v in config.ZONAL_VARS if v['name'] in ds.data_vars]
        for name in [v['name'] for v in config.ZONAL_VARS if v['name'] not in ds.data_vars]:
            print(f"  WARNING: '{name}' not in dataset, skipping zonal mean.")

        _zonal_lazy = {v['name']: compute.preload_zonal_mean(ds, v['name'])
                       for v in _zonal_valid}
        _zonal_computed = dict(zip(
            _zonal_lazy.keys(),
            dask.compute(*_zonal_lazy.values()),
        ))

        for var_cfg in _zonal_valid:
            name     = var_cfg['name']
            units    = ds[name].attrs.get('units', '')
            zonal_np = _zonal_computed[name].values   # already in memory
            for target_day in config.ZONAL_PERIODS:
                data_2d, actual_day = compute.compute_zonal_mean(
                    days, float(target_day), zonal_np)
                save_zonal_csv(data_2d, pressure_1d, lat_vals, name, actual_day)
                tag = f"{name}_day{actual_day:07.2f}"
                zonal_plots.plot_zonal_mean(
                    lat_vals, pressure_1d, data_2d, name, units,
                    actual_day, zonal_fig_dir,
                    filename=f'zonal_{tag}.png',
                    log_scale=(name in LOG_SCALE_DECADES),
                )
        _stop()
    else:
        print("  No zonal_mean_vars or zonal_mean_periods configured, skipping.")

print(f"\nDone. Figures in '{figures_dir}', data in '{data_dir}'.")

if _timing_enabled:
    _t_end_wall = _time.perf_counter()
    _t_end_dt   = _datetime.datetime.now()
    _total      = _t_end_wall - _t_start_wall
    _tracked    = sum(_timings.values())

    _right_w  = 20
    _pad      = 4   # two spaces each side
    _labels   = [f"{l}:" for l in list(_timings.keys()) + ["Other (imports/overhead)"]]
    _label_w  = max(len(l) for l in _labels + ["Started :", "Finished:", "Total   :"])
    _W        = _pad + _label_w + _right_w  # box inner width, no fixed constant

    def _row(left, right='', fill=' '):
        inner = f"  {left:<{_label_w}}{right:>{_right_w}}  "
        return f"║{inner:{fill}<{_W}}║"

    print()
    print(f"╔{'':=<{_W}}╗")
    print(f"║{'  Runtime Summary':<{_W}}║")
    print(f"╠{'':=<{_W}}╣")
    print(_row("Started :", _t_start_dt.strftime('%Y-%m-%d %H:%M:%S')))
    print(_row("Finished:", _t_end_dt.strftime('%Y-%m-%d %H:%M:%S')))
    print(_row("Total   :", f"{_total:,.1f} s"))
    print(f"╠{'':=<{_W}}╣")
    for _label, _secs in _timings.items():
        _pct = 100.0 * _secs / _total if _total > 0 else 0.0
        print(_row(f"{_label}:", f"{_secs:,.1f} s  ({_pct:5.1f}%)"))
    _other = _total - _tracked
    if _other > 0.05:
        _pct = 100.0 * _other / _total
        print(_row("Other (imports/overhead):", f"{_other:,.1f} s  ({_pct:5.1f}%)"))
    print(f"╚{'':=<{_W}}╝")

print("\n!============= Exiting exovolcano time_series diagnostics =============!")
