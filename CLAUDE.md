# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Analysis

```bash
# Run with a named experiment (looks in experiments/ directory)
python run_time_series.py ben2_vei7.yaml

# Run with an explicit path
python run_time_series.py experiments/hab2_vei7.yaml

# Run via environment variable
CONFIG=experiments/t1d_vei7.yaml python run_time_series.py

# Print usage / key reference
python run_time_series.py --help
```

Outputs go to `figures/<exp_name>/` and `data/<exp_name>/` with subdirectories `scalar/`, `profiles/`, and `aod/`. All directories are auto-created.

## Architecture

The pipeline has five modules:

**`config.py`** — Loads `experiments/<name>.yaml` at import time (via CLI arg, `CONFIG` env var, or default). Exposes path constants (`ROOT_DIR`, `FIGURES_DIR`, `DATA_DIR`), physical constants (`G_CONST`, `R_AIR`, `R_EARTH`), variable lists (`SCALAR_VARS`, `PROFILE_VARS`), and AOD parameters (`OPTICS_FILE`, `VOLC_REFF`, `RHO_AEROSOL`, `MIE_WAVE_UM`, `MIE_N_REAL`, `MIE_N_IMAG`).

**`compute.py`** — Pure computation engine; no I/O or side effects. Reads CAM NetCDF output via xarray, builds grid geometry from hybrid pressure coordinates, and computes global diagnostics. Key functions:
- `load_dataset()` — opens multi-file CAM NetCDF as a lazy xarray Dataset; extracts `gw` (Gaussian weights) from the first file separately to avoid mfdataset inflation
- `compute_geometry()` — builds pressure, layer thickness (dp), altitude (z_mid), cell area, air mass, and cell volume fields using explicit numpy broadcasting (not xarray) to avoid coordinate inflation bugs with mfdataset
- `compute_scalar()` — dispatches to `mass_integral`, `volume_integral`, or `area_mean` reduction
- `compute_profile()` — returns area-weighted mean over lat/lon, preserving `(time, lev)`

**`optics.py`** — Pure computation engine for aerosol optics; no plotting or I/O except `load_band_optics`. Key functions:
- `load_band_optics(filepath)` — opens `haze_n68_b40_mie.nc`; `rbins` are stored in cm in that file and are converted to µm on load (`* 1e4`)
- `select_band_550nm(wvn_centers)` — finds band index nearest to 18182 cm⁻¹
- `interpolate_kext(optics, i_wave, reff_um)` — log-log interpolation of Kext over rbins
- `mie_kext(wavelength_um, reff_um, ...)` — calls `miepython.efficiencies(m, diameter_um, wavelength_um)`; note this version of miepython uses diameter + wavelength in the same units, not the old size-parameter `mie(m, x)` API
- `compute_aod(volchzmd_vals, dz_m_vals, kext_cgs)` — pure numpy; multiplies dz by 100 for m→cm, sums over lev, returns `(time, lat, lon)`

**`aod_plots.py`** — AOD-specific plot functions. Callers pass the `aod/` subdirectory (`figures/<exp_name>/aod/`) so AOD figures are kept separate from top-level quicklook plots.

**`run_time_series.py`** — Orchestrator. Calls config → compute → optics → saves CSVs → makes quicklook plots. The `--help`/`-h` check happens before `config.py` is imported so it works without a YAML argument.

## Experiment YAML Schema

A fully-annotated template is at `experiments/template.yaml`. Key fields:

```yaml
root_dir: '/path/to/runs/'      # base directory for NetCDF input files
folder: 'exp_name'              # subdirectory under root_dir
file_pattern:                   # list of NetCDF filenames (CAM h1 history)
  - 'exp_name.cam.h1.0001-01-01-00000.nc'

g_const:  9.121824              # planet-specific gravity [m/s^2]
r_air:    188.965172522727      # gas constant for atmosphere [J/kg/K]
r_earth:  5797410.0             # planet radius [m]

scalar_vars:                    # global scalar time series to compute
  - name: SO2
    method: mass_integral       # or volume_integral or area_mean

profile_vars:                   # vertical profile (time, lev) to compute
  - name: T
    method: area_mean

# AOD — omit optics_file (or set to null) to skip the entire AOD section
optics_file: '/path/to/haze_n68_b40_mie.nc'
volc_reff:   1.0                # effective particle radius [µm]
rho_aerosol: 1.84               # bulk aerosol density [g/cm³]
# mie_wavelength_um:         0.55   # optional; requires miepython
# mie_refractive_index_real: 1.43
# mie_refractive_index_imag: 0.0
```

## Output Structure

```
data/<exp>/
    scalar/     <var>.csv                    — two columns: days, value
    profiles/   <var>.csv                    — comment lines with P/Z coords,
                                               then days + per-level columns
    aod/        aod_<tag>.csv               — two columns: days, global-mean AOD
                aod_zonal_<tag>.csv         — days + per-lat columns
figures/<exp>/
    quicklook_*.png                          — scalar and profile plots
    aod/        aod_<tag>_timeseries.png
                aod_<tag>_zonal_hovmoller.png
```

Profile CSVs write pressure and altitude coordinates as `# pressure_Pa:` and `# altitude_m:` comment lines before the column header (no NaN placeholders). Read in pandas with `pd.read_csv(path, comment='#')`.

## Key Domain Details

- **Physical constants are planet-specific**, not Earth-standard. Each experiment encodes its planet's gravity and radius (e.g., `g_const = 9.80665 * 0.93`). Do not substitute standard Earth values.
- **VOLCHZMD** is sulfate aerosol mass density (g/cm³); `volume_integral` multiplies by 1000 to convert to kg/m³ before integrating. AOD uses it directly in g/cm³ with a CGS Kext.
- **Gaussian weights (`gw`)**: CAM uses a Gaussian latitude grid. `gw` values sum to 2.0 and are used directly with xarray `.weighted()` for exact area-weighted means — no cosine-latitude approximation.
- **Hybrid pressure coordinates**: Pressure is computed from `hyam`, `hybm`, `hyai`, `hybi`, and `PS` using explicit numpy broadcasting, not xarray operations, to avoid a known mfdataset coordinate inflation bug.
- **Profile Hovmoller plots**: log-pressure y-axis, surface at bottom. Variables in `LOG_SCALE_DECADES` (SO2, H2SO4, Q, VOLCHZMD) use `LogNorm` colormaps; others use linear with 2nd–98th percentile clipping.
- **`rbins` unit bug in optics file**: `haze_n68_b40_mie.nc` labels `rbins` as microns but the values are in centimetres. `load_band_optics` applies `* 1e4` on load to correct this.
- **miepython API**: The installed version uses `miepython.efficiencies(m, diameter, wavelength)` with both lengths in the same units. The old `miepython.mie(m, x)` size-parameter API does not exist in this version.
