"""
compute.py - Data reading and computation engine for exovolcano analysis.

All functions are pure: no plotting, no file I/O, no side effects.

Notes on grid treatment:
  - CAM uses a Gaussian latitude grid. Gaussian weights (gw) encode the
    exact quadrature weight for each latitude band, including the poles.
    Cell area = (2*pi*R^2 / nlon) * gw(lat) is exact - no double-counting.
  - Pressure arrays are built from hybrid coefficients using explicit numpy
    broadcasting to avoid mfdataset coordinate dimension inflation bugs.
  - dp_pa is defined as p_lower_interface - p_upper_interface, which is
    always positive. All derived quantities (air mass, dz) are therefore
    also positive without any sign gymnastics.

Public API:
    load_dataset(file_list)                  -> (xr.Dataset, gw_values)
    compute_area_mean(da, gw_values)         -> DataArray (lev, or scalar per time)
    compute_geometry(ds, gw_values)          -> dict
    compute_scalar(ds, geom, name, method)   -> DataArray(time)
    compute_profile(ds, geom, name)          -> DataArray(time, lev)
    days_since_start(ds)                     -> np.ndarray(float)
"""

import numpy as np
import xarray as xr
from config import G_CONST, R_AIR, R_EARTH


# ---------------------------------------------------------------------------
# Dataset loading
# ---------------------------------------------------------------------------

def load_dataset(file_list):
    """
    Open a list of CAM NetCDF files as a single xarray Dataset.
    Loads gw cleanly from the first file to avoid mfdataset inflation.

    Returns
    -------
    ds        : xr.Dataset  (caller must close)
    gw_values : np.ndarray (nlat,), guaranteed 1D
    """
    with xr.open_dataset(file_list[0], engine='netcdf4') as ds0:
        gw_values = ds0['gw'].values.ravel()

    ds = xr.open_mfdataset(
        file_list,
        combine='by_coords',
        parallel=False,
        chunks={'time': 100},
        engine='netcdf4',
        data_vars='minimal',
    )
    print(f"Loaded {len(file_list)} file(s). Times: {len(ds['time'])}")
    return ds, gw_values


# ---------------------------------------------------------------------------
# Area-weighted global mean
# ---------------------------------------------------------------------------

def compute_area_mean(da, gw_values):
    """
    Area-weighted global mean over (lat, lon).

    Gaussian weights (gw_values) are the correct quadrature weights for
    each latitude band on the Gaussian grid. Passing them directly to
    xarray's .weighted() gives an exact area-weighted mean with no
    double-counting at the poles.

    Parameters
    ----------
    da        : DataArray with 'lat' and 'lon' dimensions
    gw_values : np.ndarray (nlat,)

    Returns
    -------
    DataArray with lat/lon reduced, other dims (time, lev) retained.
    """
    weights   = xr.DataArray(gw_values, coords={'lat': da['lat']}, dims=['lat'])
    mean_dims = [d for d in ['lat', 'lon'] if d in da.dims]
    return da.weighted(weights).mean(dim=mean_dims, skipna=True)


# ---------------------------------------------------------------------------
# Grid geometry
# ---------------------------------------------------------------------------

def compute_geometry(ds, gw_values):
    """
    Compute all grid geometry needed for mass and area integrals.

    Pressure arrays are built with explicit numpy broadcasting from hybrid
    coefficients loaded as plain 1D arrays, bypassing mfdataset inflation.

    dp_pa = p_lower_interface - p_upper_interface  (always positive)
    air_mass_cell = (dp_pa / g) * cell_area        (always positive)
    dz = R_AIR * T * dp_pa / (g * mid_p)          (always positive)

    Returns dict:
        gw_values     : np.ndarray(nlat,)
        cell_area_2d  : DataArray(lat, lon)              [m^2]
        dp_pa         : DataArray(time, lev, lat, lon)   [Pa]  positive
        mid_p         : DataArray(time, lev, lat, lon)   [Pa]
        air_mass_cell : DataArray(time, lev, lat, lon)   [kg]  positive
        dz            : DataArray(time, lev, lat, lon)   [m]   positive
        z_mid         : DataArray(time, lev, lat, lon)   [m]
        cell_volume   : DataArray(time, lev, lat, lon)   [m^3] positive
    """
    nlon = len(ds['lon'])

    # --- Cell area [m^2] ---
    # gw sums to 2.0 on [-1,1]; 2*pi*R^2 * gw / nlon gives area per cell.
    nlat = len(gw_values)
    area_band    = 2.0 * np.pi * R_EARTH**2 * gw_values        # (nlat,)
    cell_area_np = np.outer(area_band / nlon,
                            np.ones(nlon))                      # (nlat, nlon)
    cell_area_2d = xr.DataArray(cell_area_np,
                                coords={'lat': ds['lat'], 'lon': ds['lon']},
                                dims=['lat', 'lon'],
                                attrs={'units': 'm^2', 'long_name': 'grid cell area'})

    # --- Hybrid coordinate arrays ---
    # These are time-invariant 1D fields. mfdataset does not inflate them
    # since they are identical across files. .ravel() guarantees 1D numpy.
    hyam = ds['hyam'].values.ravel()   # (lev,)
    hybm = ds['hybm'].values.ravel()
    hyai = ds['hyai'].values.ravel()   # (ilev,) = (lev+1,)
    hybi = ds['hybi'].values.ravel()
    P0   = float(np.asarray(ds['P0'].values).flat[0])

    # --- Pressure fields [Pa] via explicit numpy broadcasting ---
    PS   = ds['PS'].values   # (time, lat, lon)

    # mid_p [Pa]: (time, lev, lat, lon)
    mid_p_np = (hyam[np.newaxis, :, np.newaxis, np.newaxis] * P0
                + hybm[np.newaxis, :, np.newaxis, np.newaxis] * PS[:, np.newaxis, :, :])

    # ilev_p [Pa]: (time, ilev, lat, lon)
    ilev_p_np = (hyai[np.newaxis, :, np.newaxis, np.newaxis] * P0
                 + hybi[np.newaxis, :, np.newaxis, np.newaxis] * PS[:, np.newaxis, :, :])

    # dp_pa [Pa]: lower interface minus upper interface → always positive
    # layer k: lower = ilev_p[:, k+1, :, :], upper = ilev_p[:, k, :, :]
    dp_pa_np = ilev_p_np[:, 1:, :, :] - ilev_p_np[:, :-1, :, :]   # (time, lev, lat, lon)

    coords4d = dict(time=ds['time'], lev=ds['lev'], lat=ds['lat'], lon=ds['lon'])
    dims4d   = ['time', 'lev', 'lat', 'lon']

    mid_p = xr.DataArray(mid_p_np,  coords=coords4d, dims=dims4d,
                         attrs={'units': 'Pa', 'long_name': 'layer midpoint pressure'})
    dp_pa = xr.DataArray(dp_pa_np,  coords=coords4d, dims=dims4d,
                         attrs={'units': 'Pa', 'long_name': 'layer pressure thickness (positive)'})

    # --- Air mass per cell [kg] ---
    air_mass_cell = (dp_pa / G_CONST) * cell_area_2d   # (dp_pa > 0, area > 0 → positive)

    # --- Layer thickness dz [m] from ideal gas law ---
    T_np    = ds['T'].values   # (time, lev, lat, lon)
    dz_np   = (R_AIR * T_np * dp_pa_np) / (G_CONST * mid_p_np)   # positive

    dz = xr.DataArray(dz_np, coords=coords4d, dims=dims4d,
                      attrs={'units': 'm', 'long_name': 'layer thickness'})

    # --- Midpoint altitude z_mid [m] via upward cumsum ---
    # lev index 0 = model top, nlev-1 = lowest layer (near surface).
    # Flip so index 0 = surface, cumsum upward, flip back.
    dz_flip  = dz_np[:, ::-1, :, :]                               # surface first
    z_iface  = np.concatenate([
        np.zeros((dz_np.shape[0], 1, dz_np.shape[2], dz_np.shape[3])),
        np.cumsum(dz_flip, axis=1)
    ], axis=1)                                                     # (time, nlev+1, lat, lon)
    z_mid_np = 0.5 * (z_iface[:, :-1, :, :] + z_iface[:, 1:, :, :])
    z_mid_np = z_mid_np[:, ::-1, :, :]                            # restore top-first

    z_mid = xr.DataArray(z_mid_np, coords=coords4d, dims=dims4d,
                         attrs={'units': 'm', 'long_name': 'layer midpoint altitude'})

    cell_volume = cell_area_2d * dz   # [m^3]

    # --- Diagnostics ---
    ps_mean    = float(PS[0].mean())
    z_sfc      = float(z_mid_np[0, -1, :, :].mean())
    z_top      = float(z_mid_np[0,  0, :, :].mean())
    total_area = float(cell_area_np.sum())
    expect_area = 4.0 * np.pi * R_EARTH**2

    # Mass conservation check: sum(dp_k * area / g) should equal sum(PS * area / g)
    # i.e. sum of layer masses = total atmospheric mass
    col_dp_sum  = dp_pa_np[0].sum(axis=0)              # (nlat, nlon): sum of dp over lev
    mass_from_dp = (col_dp_sum * cell_area_np / G_CONST).sum()
    mass_from_ps = (PS[0] * cell_area_np / G_CONST).sum()

    print(f"  Geometry: mean(PS)={ps_mean:.1f} Pa | z_sfc={z_sfc:.0f} m | z_top={z_top:.0f} m")
    print(f"  Area check : {total_area:.6e} m^2  (expect {expect_area:.6e},"
          f" ratio={total_area/expect_area:.6f})")
    print(f"  Mass check : from dp={mass_from_dp:.6e} kg,"
          f" from PS={mass_from_ps:.6e} kg,"
          f" ratio={mass_from_dp/mass_from_ps:.6f}")

    return {
        'gw_values':     gw_values,
        'cell_area_2d':  cell_area_2d,
        'dp_pa':         dp_pa,
        'mid_p':         mid_p,
        'air_mass_cell': air_mass_cell,
        'dz':            dz,
        'z_mid':         z_mid,
        'cell_volume':   cell_volume,
    }


# ---------------------------------------------------------------------------
# Scalar time series
# ---------------------------------------------------------------------------

def compute_scalar(ds, geom, name, method):
    """
    Compute a global scalar time series for one variable.

    mass_integral   : total_mass [kg] = sum_over_(lev,lat,lon)[ q * air_mass_cell ]
    volume_integral : total_mass [kg] = sum_over_(lev,lat,lon)[ rho * cell_volume ]
                      (VOLCHZMD is g/cm^3, converted to kg/m^3 * m^3 = kg)
    area_mean       : area-weighted global mean via compute_area_mean()

    Parameters
    ----------
    ds     : xr.Dataset
    geom   : dict from compute_geometry()
    name   : str
    method : 'mass_integral' | 'volume_integral' | 'area_mean'

    Returns
    -------
    DataArray(time) or None
    """
    if name not in ds.data_vars:
        print(f"  WARNING: '{name}' not found in dataset, skipping.")
        return None

    var = ds[name]

    if method == 'mass_integral':
        ref = geom['air_mass_cell'].reindex_like(var)
        return (var * ref).sum(dim=['lev', 'lat', 'lon'], skipna=True).load()

    elif method == 'volume_integral':
        ref = geom['cell_volume'].reindex_like(var)
        scale = 1000.0 if name == 'VOLCHZMD' else 1.0   # g/cm^3 -> kg/m^3
        return (var * ref * scale).sum(dim=['lev', 'lat', 'lon'], skipna=True).load()

    elif method == 'area_mean':
        return compute_area_mean(var, geom['gw_values']).load()

    else:
        print(f"  WARNING: Unknown method '{method}' for '{name}', skipping.")
        return None


# ---------------------------------------------------------------------------
# Profile time series
# ---------------------------------------------------------------------------

def compute_profile(ds, geom, name):
    """
    Global area-weighted mean vertical profile time series.
    Returns DataArray(time, lev) or None.
    """
    if name not in ds.data_vars:
        print(f"  WARNING: '{name}' not found in dataset, skipping.")
        return None
    var = ds[name]
    if 'lev' not in var.dims:
        print(f"  WARNING: '{name}' has no 'lev' dimension, skipping profile.")
        return None
    return compute_area_mean(var, geom['gw_values']).load()


# ---------------------------------------------------------------------------
# Time coordinate utility
# ---------------------------------------------------------------------------

def days_since_start(ds):
    """
    Elapsed time in days from the first timestep.
    Works for cftime and numpy datetime64.
    """
    times = ds['time'].values
    t0    = times[0]
    try:
        return np.array([(t - t0).days + (t - t0).seconds / 86400.0 for t in times])
    except AttributeError:
        return (times - t0).astype('timedelta64[s]').astype(float) / 86400.0
