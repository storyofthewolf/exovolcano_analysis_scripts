# exovolcano_analysis_scripts

Post-processing diagnostics for CAM volcanic eruption simulations. Reads CAM `h1` NetCDF history files, computes global scalar and vertical profile time series, and produces quicklook figures. Designed for both Earth-like and exoplanet atmospheres — physical constants (gravity, gas constant, planet radius) are per-experiment rather than hard-coded.

## Requirements

- Python 3.9+
- `numpy`, `xarray`, `matplotlib`, `pyyaml`, `netcdf4`
- `miepython` — only needed if using the optional Mie AOD calculation

## Running

```bash
python run_time_series.py <experiment>.yaml
```

The script looks for `<experiment>.yaml` in the `experiments/` directory. You can also pass a full path or set the `CONFIG` environment variable:

```bash
python run_time_series.py experiments/t1d_vei7.yaml
CONFIG=experiments/ben2_vei7.yaml python run_time_series.py
```

Outputs are written under `figures/<exp_name>/` and `data/<exp_name>/`.

### Runtime flags

Any combination of these flags can be appended to the command:

| Flag | Effect |
|------|--------|
| `--time` | Print a per-section timing summary at the end |
| `--no-scalars` | Skip scalar time series (computation + CSV) |
| `--no-profiles` | Skip profile time series (computation + CSV) |
| `--no-plots` | Skip all figure output (CSV data is still written) |
| `--no-aod` | Skip AOD calculation |
| `--no-zonal` | Skip zonal mean snapshots |

These flags override the YAML config — you do not need to comment out YAML sections to skip a step. They are particularly useful on the cluster where zonal mean snapshots are expensive:

```bash
# Skip the two most expensive sections for a fast diagnostic pass
python run_time_series.py ben2_vei7.yaml --no-zonal --no-aod --time

# Write CSVs only, no figures
python run_time_series.py ben2_vei7.yaml --no-plots
```

## Experiments

**All configuration — including AOD settings — lives in the experiment YAML file** passed on the command line. There is no separate config file to edit. To enable a feature, add the corresponding keys to your YAML.

A fully-annotated template covering every available field is at:

```
experiments/template.yaml
```

Copy it, rename it for your run, and fill in your values.

Existing experiment configs:

| File | Planet / run |
|------|-------------|
| `t1d_vei7.yaml` | TRAPPIST-1d, VEI-7 eruption |
| `hab1_vei7.yaml` | Habitable zone planet 1, VEI-7 |
| `hab2_vei7.yaml` | Habitable zone planet 2, VEI-7 |
| `ben2_vei7.yaml` | Benchmark planet 2, VEI-7 (exo-atmosphere constants) |
| `modern_vei6.yaml` | Modern Earth, VEI-6 (terrestrial constants) |

## Experiment YAML reference

The sections below describe every key. All optional keys have defaults and can be omitted.

### Input files (required)

```yaml
root_dir: '/path/to/cesm/archive/'
folder:   'run_name/atm/hist'        # subdirectory under root_dir
file_pattern:
  - 'run_name.cam.h1.0001-01-01-00000.nc'
  - 'run_name.cam.h1.0001-04-11-00000.nc'
  # one line per CAM h1 file to include
```

The filename prefix before the first `.` is used as the experiment name for output subdirectory naming, so all files must share the same prefix.

### Physical constants (required)

These must match the planet being simulated. Do not use Earth values for exoplanet runs — grid geometry, layer mass, and AOD all depend on them.

```yaml
g_const: 9.80665         # surface gravity [m/s²]  Earth: 9.80665
r_air:   287.058         # specific gas constant [J/kg/K]  Earth air: 287, pure CO2: 188.965
r_earth: 6.371e6         # mean planet radius [m]
```

### Output directories (optional)

```yaml
figures_dir: 'figures'   # PNGs go to figures/<exp_name>/
data_dir:    'data'      # CSVs go to data/<exp_name>/scalar/ and .../profiles/
```

### Scalar time series (optional)

```yaml
scalar_vars:
  - name: SO2
    method: mass_integral     # Σ(q × air_mass_cell) [kg]
  - name: VOLCHZMD
    method: volume_integral   # Σ(ρ × cell_volume) [kg]; VOLCHZMD g/cm³→kg/m³ handled internally
  - name: TS
    method: area_mean         # Gaussian-weighted global mean [native units]
```

### Vertical profile time series (optional)

```yaml
profile_vars:
  - name: T
    method: area_mean    # always area_mean regardless of value given
```

### Zonal mean snapshots (optional)

Produces a contour plot and CSV for each variable at each requested day. The nearest available timestep is used.

```yaml
zonal_mean_vars:
  - name: T
  - name: SO2
  - name: H2SO4
  - name: VOLCHZMD

zonal_mean_periods:
  increment: [0, 1, 4, 10, 50, 100]   # days since start
```

Both sections must be present to enable zonal output. Omitting either silently skips the section.

### AOD calculation (optional)

AOD requires VOLCHZMD in the CAM output. **All AOD settings are in the same experiment YAML.** If `optics_file` is absent or `null`, the entire AOD section is skipped with a printed warning.

```yaml
# Required to enable AOD:
optics_file: '/path/to/haze_n68_b40_mie.nc'   # pre-computed Mie optics table

# Optional tuning (defaults shown):
volc_reff:   1.0    # effective particle radius [µm] for Kext lookup/Mie
rho_aerosol: 1.84   # bulk aerosol density [g/cm³] — used only by Mie path

# Optional single-wavelength Mie AOD (uncomment all three to enable):
# mie_wavelength_um:         0.55   # wavelength [µm]
# mie_refractive_index_real: 1.43   # real refractive index
# mie_refractive_index_imag: 0.0    # imaginary part (positive convention)
```

The band-interpolated 550 nm AOD is always computed when `optics_file` is set. The Mie block adds a second AOD calculation at a custom wavelength and requires `miepython` to be installed.

---

## Plots produced

### Scalar time series

| Filename | Contents |
|----------|----------|
| `quicklook_sulfur_budget.png` | SO2 + H2SO4 + VOLCHZMD total masses plus their sum |
| `quicklook_TS.png` | Global mean surface temperature |
| `quicklook_TGCLDLWP.png` | Global mean total liquid water path |
| `quicklook_TMQ.png` | Global mean total precipitable water |
| `quicklook_Q_mass.png` | Global total water vapor mass |

Each plot uses a shared x-axis of days since the first time step.

### Profile Hovmoller diagrams

One `quicklook_profile_<VAR>.png` per entry in `profile_vars`. The y-axis is log-pressure in mb with the surface at the bottom. Variables in the log-scale group (SO2, H2SO4, Q, VOLCHZMD) use a `LogNorm` colormap spanning 4 orders of magnitude below the peak value. All other variables use a linear colormap clipped to the 2nd–98th percentile of the data.

### Zonal mean contour plots (requires `zonal_mean_vars` and `zonal_mean_periods`)

One `zonal/zonal_<VAR>_day<DAY>.png` per (variable, day) pair. The x-axis is latitude in degrees (full range), the y-axis is log-pressure in mb with the surface at the bottom — the same layout as the profile Hovmoller plots. Log-scale variables (SO2, H2SO4, Q, VOLCHZMD) use `LogNorm` spanning 4 decades below the peak; all others use linear scale clipped to the 2nd–98th percentile.

### AOD plots (requires `optics_file`)

Two plots are produced per wavelength:

| Filename | Contents |
|----------|----------|
| `quicklook_aod_550nm_band_timeseries.png` | Global mean AOD at 550 nm (band-interpolated Kext) |
| `quicklook_aod_550nm_band_zonal_hovmoller.png` | Zonal mean AOD vs time (linear latitude axis) |
| `quicklook_aod_<W>um_mie_timeseries.png` | Global mean AOD at `mie_wavelength_um` (Mie Kext) |
| `quicklook_aod_<W>um_mie_zonal_hovmoller.png` | Zonal mean AOD at `mie_wavelength_um` vs time |

The Mie pair is only produced when `mie_wavelength_um` is set in the YAML.

---

## Zonal mean CSV format

Files are written to `data/<exp>/zonal/<VAR>_day<DAY>.csv`. Rows are pressure levels, columns are latitudes:

```
pressure_mb, -90.0000, -87.5000, ..., 90.0000
1.2345,       0.0,      0.0,     ..., 0.0
2.5678,       0.0,      0.0,     ..., 0.0
...
```

Read in pandas with `pd.read_csv(path)`. The `pressure_mb` column gives the layer midpoint pressure in mb (time- and area-mean from the first timestep).

---

## AOD calculation details

AOD is computed from the VOLCHZMD field (aerosol mass density, g/cm³) and the layer thickness dz (m):

```
AOD_layer = Kext [cm²/g] × ρ [g/cm³] × dz [cm]
AOD_column = Σ AOD_layer  (sum over lev)
```

### Band-interpolated Kext (default)

The optics file `haze_n68_b40_mie.nc` contains pre-computed Mie extinction efficiencies on a grid of 68 wavenumber bands × 40 particle radii. The 550 nm band is found by locating the band center nearest to 18182 cm⁻¹. Kext is then interpolated to `volc_reff` using log-log interpolation over the radius axis.

### Single-wavelength Mie Kext (optional)

When `mie_wavelength_um` is set, `miepython` computes Kext directly:

```
x      = 2π × r_eff / λ
m      = n_real - i × |n_imag|    (miepython sign convention)
Kext   = Q_ext × 3 / (4 × r_eff_cm × ρ_bulk)   [cm²/g]
```

This is useful for wavelengths outside the optics table or for sensitivity tests with different refractive indices.

---

## Running on an HPC cluster

### Login node vs. batch job

**Avoid running on the login node for full production runs.** Login nodes are shared across all users on the cluster. Running a compute-heavy script there means:

- Your dask threads compete with other users' processes for CPU time, causing 2–4× run-to-run timing variance even with identical code and data.
- Sustained CPU usage on a login node is against policy on most HPC systems (SLURM, PBS, LSF) and may result in your process being killed or your account flagged.
- Lustre I/O from the login node is uncontested for bandwidth but the metadata server is shared — opening many files in parallel can degrade performance for other users.

Use the login node only for quick tests (`--no-zonal --no-aod`) and for developing/debugging. Submit a batch job for full runs.

### Controlling dask thread count

The script uses a threaded dask scheduler. The number of workers defaults to `min(8, cpu_count)` but can be overridden with the `DASK_WORKERS` environment variable:

```bash
# Login node — keep low to avoid contention
DASK_WORKERS=4 python run_time_series.py ben2_vei7.yaml --time

# Batch job — set to your allocated core count
DASK_WORKERS=32 python run_time_series.py ben2_vei7.yaml --time
```

### Example SLURM batch script

Save as `submit_analysis.sh` and submit with `sbatch submit_analysis.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=exovol_analysis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # dask threads — adjust to what you need
#SBATCH --mem=64G                # set based on dataset size; 32–64 GB typical
#SBATCH --time=01:00:00          # walltime limit
#SBATCH --output=logs/%j_analysis.out
#SBATCH --error=logs/%j_analysis.err

# Activate your Python environment
source activate exovol            # or: module load python; source /path/to/venv/bin/activate

# Tell the script how many cores it has
export DASK_WORKERS=$SLURM_CPUS_PER_TASK

# Run the analysis
python run_time_series.py ben2_vei7.yaml --time
```

Adjust `--cpus-per-task` and `--mem` to match your dataset. A 1000-timestep run at 96×144 grid with 5 scalar variables and 5 profile variables typically needs 32–64 GB and runs in under 60 s with 16 dedicated cores.

### Example PBS/Torque batch script

```bash
#!/bin/bash
#PBS -N exovol_analysis
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=01:00:00
#PBS -o logs/analysis.out
#PBS -e logs/analysis.err

cd $PBS_O_WORKDIR
source activate exovol

export DASK_WORKERS=16
python run_time_series.py ben2_vei7.yaml --time
```

### Quick interactive job (srun)

For interactive testing with a dedicated core allocation:

```bash
srun --ntasks=1 --cpus-per-task=8 --mem=32G --time=00:30:00 --pty bash
source activate exovol
export DASK_WORKERS=8
python run_time_series.py ben2_vei7.yaml --no-zonal --time
```

---

## Code structure

```
run_time_series.py   Orchestrator: loads data, calls compute/optics, saves outputs
config.py            YAML loader; exposes all experiment parameters as module constants
compute.py           Pure computation: grid geometry, mass integrals, area means, zonal means
optics.py            Pure computation: optics table I/O, Kext interpolation, Mie, AOD
aod_plots.py         AOD-specific plot functions (timeseries, zonal Hovmoller)
zonal_plots.py       Zonal mean contour plot function; defines LOG_SCALE_DECADES
experiments/         One YAML per model run
data/                CSV output (scalar/, profiles/, aod/, zonal/ subdirectories per experiment)
figures/             PNG output (one subdirectory per experiment, with aod/ and zonal/ sub-dirs)
```

`compute.py` and `optics.py` contain no plotting or file I/O (except `optics.load_band_optics`). All grid geometry — pressure from hybrid coefficients, layer thickness, cell area, air mass — is computed in `compute.compute_geometry()` and passed downstream as a plain dict of DataArrays.
