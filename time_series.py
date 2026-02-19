import xarray as xr
import matplotlib.pyplot as plt
import os
import numpy as np
import glob
import config
from config import G_CONST, R_AIR, R_EARTH, OUTPUT_DIR

file_list = config.get_file_list()
print(f"file_list:  '{file_list}'")

experiment_name = config.get_experiment_name()


print(f"experiment name:  '{experiment_name}'")

def get_units(data_array):
    """Safely retrieves units attribute from a DataArray."""
    return data_array.attrs.get('units', 'unknown units')

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

try:
    # Use open_mfdataset for multi-file support
    # Added 'data_vars="minimal"' to reduce metadata overhead
    with xr.open_mfdataset(
        file_list, 
        combine='by_coords', 
        parallel=False,
        chunks={'time': 100},
        engine='netcdf4',
        data_vars='minimal'
    ) as ds:
        
        print("All datasets loaded successfully.")

        # --- NaN Diagnostics ---
        print("\n--- Running NaN Diagnostics ---")
        for var in ['PS', 'T', 'SO2']:
            if var in ds.data_vars:
                nan_count = ds[var].isnull().sum().compute()
                print(f"Checking '{var}': Found {nan_count} NaN values.")
        print("--- End NaN Diagnostics ---\n")
        
        # --- 1. Calculate Grid Geometry ---
        print("Calculating grid geometry...")
        # Gaussian weights (gw) represent the area of each latitude band
        with xr.open_dataset(file_list[0], engine='netcdf4') as first_ds:
            gw_values = first_ds['gw'].load().values
            
        area_band_values = 2.0 * np.pi * (R_EARTH**2) * gw_values
        nlon = len(ds['lon'])
        cell_area_1d_da = xr.DataArray(area_band_values / nlon, coords={'lat': ds['lat']}, dims=['lat'])
        ones_2d = xr.ones_like(ds['PS'].isel(time=0, drop=True))
        cell_area_2d = (cell_area_1d_da * ones_2d)

        # Pressure calculations for hybrid levels
        mid_pressure_pa = ds['hyam'] * ds['P0'] + ds['hybm'] * ds['PS']
        ilev_pressure_pa = ds['hyai'] * ds['P0'] + ds['hybi'] * ds['PS']
        dp_pa = ilev_pressure_pa.diff(dim='ilev').rename({'ilev': 'lev'})
        dp_pa['lev'] = ds['lev']

        # --- 2. Calculate Masses ---
        air_mass_cell = (dp_pa / G_CONST) * cell_area_2d
        
        cell_volume = None
        if 'VOLCHZMD' in ds.data_vars and 'T' in ds.data_vars:
            cell_height_dz = (dp_pa * R_AIR * ds['T']) / (mid_pressure_pa * G_CONST)
            cell_volume = cell_area_2d * cell_height_dz

        # --- 3. Compute Time Series ---
        print("Calculating total mass for each variable...")
        variables_to_plot = {}
        
        def compute_mass_ts(var_name, mass_ref):
            if var_name in ds.data_vars:
                print(f"Processing {var_name}...")
                aligned_ref = mass_ref.reindex_like(ds[var_name])
                return (ds[var_name] * aligned_ref).sum(dim=['lev', 'lat', 'lon'], skipna=True).load()
            return None

        variables_to_plot['SO2'] = compute_mass_ts('SO2', air_mass_cell)
        variables_to_plot['H2SO4'] = compute_mass_ts('H2SO4', air_mass_cell)
        variables_to_plot['Q'] = compute_mass_ts('Q', air_mass_cell)
        
        if 'VOLCHZMD' in ds.data_vars and cell_volume is not None:
            # Convert VOLCHZMD from g/cm3 to kg/m3 (*1000)
            variables_to_plot['VOLCHZMD'] = compute_mass_ts('VOLCHZMD', cell_volume * 1000.0)

        # Filter missing
        variables_to_plot = {k: v for k, v in variables_to_plot.items() if v is not None}

        # Sum Sulfur components
        s_vars = ['SO2', 'H2SO4', 'VOLCHZMD']
        available_s = [variables_to_plot[v] for v in s_vars if v in variables_to_plot]
        if available_s:
            variables_to_plot['Total Sulfur (SO2+H2SO4+VOLCHZMD)'] = sum(available_s)

        # --- 4. Plotting Logic ---
        print("Generating plots...")
        
        def prep_for_plot(da):
            try:
                t0 = ds['time'].values[0]
                days = np.array([(t - t0).days + (t - t0).seconds/86400.0 for t in ds['time'].values])
                da_plot = da.copy()
                da_plot['time'] = days
                return da_plot, "Days since start"
            except Exception:
                return da, ds['time'].attrs.get('units', 'time')

        if variables_to_plot:
            plot_order = ['SO2', 'H2SO4', 'VOLCHZMD', 'Total Sulfur (SO2+H2SO4+VOLCHZMD)', 'Q']
            ordered_keys = [k for k in plot_order if k in variables_to_plot]
            
            fig, axes = plt.subplots(nrows=len(ordered_keys), ncols=1, figsize=(12, len(ordered_keys)*3), sharex=True)
            if len(ordered_keys) == 1: axes = [axes]
            
            for i, key in enumerate(ordered_keys):
                data_to_plot, time_label = prep_for_plot(variables_to_plot[key])
                data_to_plot.plot(ax=axes[i])
                axes[i].set_title(f"Global Total Mass of {key}")
                axes[i].set_ylabel("Mass (kg)")
                axes[i].grid(True)
                axes[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            axes[-1].set_xlabel(time_label)
            plt.tight_layout()
            filename = f"{experiment_name}.global_mean_timeseries_v1.png"
            save_fig1 = os.path.join(OUTPUT_DIR, filename)
            plt.savefig(save_fig1)
            plt.close()

        # --- 5. Global Mean TS ---
        mean_vars = ['TS', 'TGCLDLWP', 'TMQ']
        mean_data = {}
        weights_da = xr.DataArray(gw_values, coords={'lat': ds['lat']}, dims=['lat'])
        
        for v in mean_vars:
            if v in ds.data_vars:
                print(f"Calculating mean for {v}...")
                mean_data[v] = ds[v].weighted(weights_da).mean(dim=('lat', 'lon'), skipna=True).load()

        if mean_data:
            fig, axes = plt.subplots(nrows=len(mean_data), ncols=1, figsize=(12, len(mean_data)*3), sharex=True)
            if len(mean_data) == 1: axes = [axes]
            
            for i, (key, val) in enumerate(mean_data.items()):
                data_to_plot, time_label = prep_for_plot(val)
                data_to_plot.plot(ax=axes[i])
                axes[i].set_title(f"Global Mean {key}")
                axes[i].set_ylabel(get_units(ds[key]))
                axes[i].grid(True)
            
            axes[-1].set_xlabel(time_label)
            plt.tight_layout()
            filename = f"{experiment_name}.global_mean_timeseries_v2.png"
            save_fig2 = os.path.join(OUTPUT_DIR, filename)
            plt.savefig(save_fig2)
            plt.close()

except Exception as e:
    print(f"An error occurred: {e}")
    import traceback
    traceback.print_exc()

print(f"\nScript finished. Plots are in '{OUTPUT_DIR}'.")
