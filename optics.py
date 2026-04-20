"""
optics.py - Aerosol optical property functions for exovolcano analysis.

Pure functions only. No plotting, no I/O except load_band_optics.

Public API:
    load_band_optics(filepath)              -> dict of numpy arrays
    select_band_550nm(wvn_centers)          -> int
    interpolate_kext(optics, i_wave, reff)  -> float [cm²/g]
    mie_kext(wavelength_um, reff_um, ...)   -> float [cm²/g]
    compute_aod(volchzmd_vals, dz_m_vals, kext_cgs) -> np.ndarray (time, lat, lon)
"""

import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Band optics file loader
# ---------------------------------------------------------------------------

def load_band_optics(filepath):
    """
    Open haze_n68_b40_mie.nc and return a dict of numpy arrays.

    Returns
    -------
    dict with keys:
        wvnrng      : (69,)    wavenumber edges [cm⁻¹]
        wvn_centers : (68,)    midpoints of adjacent wvnrng edges [cm⁻¹]
        rbins       : (40,)    particle radii from rbins[:, 0] [microns]
        Kext        : (68, 40) extinction efficiency [cm²/g]
        W           : (68, 40) single scattering albedo
        G           : (68, 40) asymmetry parameter
    """
    with xr.open_dataset(filepath, engine='netcdf4') as ds:
        wvnrng      = ds['wvnrng'].values.ravel()           # (69,)
        wvn_centers = 0.5 * (wvnrng[:-1] + wvnrng[1:])     # (68,)
        rbins       = ds['rbins'].values[:, 0] * 1e4           # cm -> µm
        Kext        = ds['Kext'].values[:, :, 0]             # (68, 40)
        W           = ds['W'].values[:, :, 0]                # (68, 40)
        G           = ds['G'].values[:, :, 0]                # (68, 40)

    return {
        'wvnrng':      wvnrng,
        'wvn_centers': wvn_centers,
        'rbins':       rbins,
        'Kext':        Kext,
        'W':           W,
        'G':           G,
    }


# ---------------------------------------------------------------------------
# Band selection
# ---------------------------------------------------------------------------

def select_band_550nm(wvn_centers):
    """
    Return the integer index of the band whose center is closest to
    18182 cm⁻¹ (= 1 / 0.55e-4 cm = 550 nm).

    Parameters
    ----------
    wvn_centers : array-like (nbands,) [cm⁻¹]

    Returns
    -------
    int
    """
    target = 1.0 / 0.55e-4   # 18181.8... cm⁻¹
    return int(np.argmin(np.abs(np.asarray(wvn_centers) - target)))


# ---------------------------------------------------------------------------
# Kext interpolation from band optics table
# ---------------------------------------------------------------------------

def interpolate_kext(optics, i_wave, reff_um):
    """
    Interpolate Kext[i_wave, :] over rbins at reff_um using log-log
    interpolation.

    Parameters
    ----------
    optics   : dict from load_band_optics()
    i_wave   : int   wavelength band index
    reff_um  : float effective radius [microns]

    Returns
    -------
    float  Kext [cm²/g]

    Raises
    ------
    ValueError if reff_um is outside the range of rbins.
    """
    rbins = optics['rbins']
    kext  = optics['Kext'][i_wave, :]

    r_min, r_max = rbins.min(), rbins.max()
    if reff_um < r_min or reff_um > r_max:
        raise ValueError(
            f"reff_um={reff_um} is outside the optics table range "
            f"[{r_min}, {r_max}] microns."
        )

    log_r    = np.log(rbins)
    log_kext = np.log(kext)
    return float(np.exp(np.interp(np.log(reff_um), log_r, log_kext)))


# ---------------------------------------------------------------------------
# Single-wavelength Mie calculation
# ---------------------------------------------------------------------------

def mie_kext(wavelength_um, reff_um, n_real, n_imag, rho_bulk_gcc):
    """
    Compute Kext [cm²/g] for a monodisperse sphere at a single wavelength
    using Mie theory (miepython).

    Parameters
    ----------
    wavelength_um  : float  wavelength [microns]
    reff_um        : float  effective radius [microns]
    n_real         : float  real part of refractive index
    n_imag         : float  imaginary part (positive convention; sign is
                            handled internally per miepython convention)
    rho_bulk_gcc   : float  bulk aerosol density [g/cm³]

    Returns
    -------
    float  Kext [cm²/g]
    """
    import miepython

    # miepython.efficiencies takes diameter and wavelength in the same units;
    # pass both in µm.  Negative imaginary part is the miepython convention.
    m = complex(n_real, -abs(n_imag))
    diameter_um = 2.0 * reff_um
    Qext, _Qsca, _Qback, _g = miepython.efficiencies(m, diameter_um, wavelength_um)

    reff_cm = reff_um * 1e-4
    kext_cgs = Qext * 3.0 / (4.0 * reff_cm * rho_bulk_gcc)
    return float(kext_cgs)


# ---------------------------------------------------------------------------
# AOD computation
# ---------------------------------------------------------------------------

def compute_aod(volchzmd_vals, dz_m_vals, kext_cgs):
    """
    Compute column aerosol optical depth from 3D aerosol density and layer
    thickness.

    AOD per layer = kext_cgs * volchzmd_vals * (dz_m_vals * 100)
    where the factor of 100 converts dz from metres to centimetres.
    VOLCHZMD is already in g/cm³, so no unit conversion is needed for it.

    Parameters
    ----------
    volchzmd_vals : np.ndarray (time, lev, lat, lon)  [g/cm³]
    dz_m_vals     : np.ndarray (time, lev, lat, lon)  [m]
    kext_cgs      : float                             [cm²/g]

    Returns
    -------
    np.ndarray (time, lat, lon)  column AOD [dimensionless]
    """
    dz_cm   = dz_m_vals * 100.0                                 # m -> cm
    aod_3d  = kext_cgs * volchzmd_vals * dz_cm                  # (time, lev, lat, lon)
    return aod_3d.sum(axis=1)                                    # sum over lev
