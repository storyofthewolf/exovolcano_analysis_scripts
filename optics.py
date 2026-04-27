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
    Open volc_pw1975_n68_r1.0um_mie.nc and return a dict of numpy arrays.

    NOTE: Currently only supports optics files with nbins=1 (single particle
    radius). Multi-bin files will load without error but interpolation over
    rbins has not been validated — additional work would be needed to support
    multiple radii.

    Returns
    -------
    dict with keys:
        wvnrng      : (69,)    wavenumber edges [cm⁻¹]
        wvn_centers : (68,)    midpoints of adjacent wvnrng edges [cm⁻¹]
        rbins       : (nbins,) particle radii [cm] as stored in the file
        Kext        : (68, nbins) extinction efficiency [cm²/g]
        W           : (68, nbins) single scattering albedo
        G           : (68, nbins) asymmetry parameter
    """
    with xr.open_dataset(filepath, engine='netcdf4') as ds:
        wvnrng      = ds['wvnrng'].values.ravel()           # (69,)
        wvn_centers = 0.5 * (wvnrng[:-1] + wvnrng[1:])     # (68,)
        rbins       = ds['rbins'].values[:, 0]              # (nbins,) in cm
        kext_raw    = ds['Kext'].values[:, :, 0]            # (68, nbins)
        kext_units  = ds['Kext'].attrs.get('units', '').strip()
        W           = ds['W'].values[:, :, 0]               # (68, nbins)
        G           = ds['G'].values[:, :, 0]               # (68, nbins)

    if kext_units in ('m2 kg-1', 'm2/kg', 'm^2/kg', 'm2 kg^-1'):
        print("-- SI units detected in optpropt file --")
        print("-- convert to cgs for calculation (x10) --")
        kext_cgs = kext_raw * 10.0                          # m²/kg -> cm²/g
    else:
        print(" CGS units detected in optpropt file ")
        kext_cgs = kext_raw                                 # already cm²/g

    return {
        'wvnrng':      wvnrng,
        'wvn_centers': wvn_centers,
        'rbins':       rbins,       # cm (as stored in file)
        'Kext':        kext_cgs,    # cm²/g
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
    Return Kext [cm²/g] for wavelength band i_wave.

    NOTE: Only nbins=1 optics files are currently supported. For single-bin
    files, reff_um is not used in the lookup — Kext is read directly from the
    one available bin. Multi-bin interpolation (log-log over rbins) is
    implemented below but has not been validated; supporting it properly would
    require confirming rbins units and testing against known references.

    Parameters
    ----------
    optics   : dict from load_band_optics()
    i_wave   : int   wavelength band index
    reff_um  : float effective radius [microns] (unused for nbins=1)

    Returns
    -------
    float  Kext [cm²/g]

    Raises
    ------
    ValueError if reff_um is outside the range of rbins (multi-bin only).
    """
    rbins = optics['rbins']
    kext  = optics['Kext'][i_wave, :]

    if len(rbins) == 1:
        # Single-bin file: radius is not a free parameter; return directly.
        return float(kext[0])

    r_min, r_max = rbins.min(), rbins.max()
    rtol = 0.01
    if reff_um < r_min * (1 - rtol) or reff_um > r_max * (1 + rtol):
        raise ValueError(
            f"reff_um={reff_um} is outside the optics table range "
            f"[{r_min:.4f}, {r_max:.4f}] (file units)."
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
