# WCB Trajectory Analysis

This repository contains scripts and tools for the analysis of Warm Conveyor Belt (WCB) trajectory data.

Created by **Cornelis Schwenk** (<corny.schwenk@gmail.com>) during work at the *Johannes Gutenberg University Mainz* in the group of **Annette Miltenberger**.

---

## 📁 File Structure

### `NAMELIST_X.jl` (X = 1, 2, 3, ...)
Defines:
- Paths to trajectory data
- Variable names
- Simulation start/end time
- Number of timesteps
- Time step interval (seconds)
- Longitude/latitude of WCB at selection times

### `CASE.jl`
- Defines the case number **X** corresponding to the selected WCB case.
- Each case must have its own path pointing to `NAMELIST_X.jl`.

### `analysis_mogon.slurm`
SLURM job script that:
- Merges trajectory files (requires **CDO - Climate Data Operators**)
- Runs `run_WCB_selection.jl`
- Uses functions from `WCB_selection.jl` to produce NetCDF files containing trajectory and microphysical variables (e.g., Condensation Ratio (CR), Precipitation Efficiency (PE))

### `run_WCB_selection.jl`
Runs the WCB selection and NetCDF creation using functions from `WCB_selection.jl`.

### `WCB_selection.jl`
Contains functions for writing NetCDF files with WCB trajectory data.

### `WCB_functions.jl`
Provides functions for trajectory data calculations, including:
- Precipitation Efficiency (PE)
- Condensation Ratio (CR)
- Hydrometeor flux
- Other microphysical diagnostics

### `use_config.jl`
Include this when using the Julia REPL for improved usability.

### `statistical_functions.jl`
Functions for statistical data processing and manipulation.

### `thermodynamic_functions.jl`
Implements thermodynamic utilities (e.g., RH\_i, RH\_w).

### `install_packages.jl`
Installs all required Julia packages.

---

## 🚀 Instructions

### 1️⃣ Install Required Packages
```bash
julia --project=. install_packages.jl
```
### 2️⃣ Specify the Case to Analyze
Edit the file `CASE.jl` and set the desired case number X.

### 3️⃣ Configure Paths and Variables
Edit `NAMELIST_X.jl` corresponding to your selected case.

### 4️⃣ Run the Selection Script
Modify paths in `analysis_mogon.slurm` and submit to SLURM
Alternatively, run the same commands manually in the terminal.

----------
## ☁️  How the WCB Selection Algorithm Works
Warm Conveyor Belt (WCB) trajectories are identified in two stages using the information specified in each `NAMELIST_X.jl` file and the functions in `WCB_selection.jl`.
The selection is based on physical ascent criteria and cyclone location checks.

### 1. Reading Trajectory and Configuration Metadata.
`WCB_selection.jl` reads the NAMELIST via:
``` 
    include("CASE.jl")
    include("NAMELIST_$(CASE).jl")
```
From the NAMELIST, it receives:
    **File paths** (`TRAJ_PATH`, `NWP_PATH`, …)
    **Time settings** (`START_OF_SIM`, `START_OF_TRAJS`, `TRAJ_TIME_STEP`)
    **Physical units** (e.g., pressure in Pa, lat/lon in radians)
    **Two cyclone "boxes"** defined by:
```
    TIME_CHECK_1, LON_CHECK_1, LAT_CHECK_1
    TIME_CHECK_2, LON_CHECK_2, LAT_CHECK_2
```
These limits are defined by you and indicate where the cyclone is located at two key times.
Look at `NAMELIST_1.jl` line 133 to see an example of how the lon and lat values create a box around the WCB. 
You have to create plots of the NWP data (cloud cover) to select these times and lon-lat values.

### 2.  Stage 1: Identifying WCB Ascent Trajectories
The first filtering uses only trajectory vertical motion.
A trajectory is marked as a WCB candidate if it **rises at least 600 hPa within ≤ 48 hours**.
This is implemented in:
```
    find_WCB(ds)
```
The key function here is: 
```
    find_tau_and_asc(p)
```
This finds:
    - `a`: when the ascent begins (index in time)
    - `tau600`: how long the 600 hPa rise takes (in time steps)

### 3. Stage 2: Cyclone Location Check
Not all ascending air parcels belong to the cyclone of interest.
To ensure the selected trajectories correspond to the desired WCB, they must be located inside two **longitude/latitude boxes** at two different times during the development of the cyclone.

Coordinates from the trajectory file are converted from radians to degrees. Then the algorithm checks whether each candidate lies within the boxes at:
```
    j_check1 = dt_to_nearest_rt(TIME_CHECK_1)
    j_check2 = dt_to_nearest_rt(TIME_CHECK_2)
```
For each candidate trajectory:
```
    if NOT in box 1 at TIME_CHECK_1 → remove
    if NOT in box 2 at TIME_CHECK_2 → remove
```

This filtering happens inside:
```
    WCB_sel(ds)
```
At the end, only WCB trajectories that both ascend rapidly AND are inside the cyclone at both times are kept.

### 4. Output of the Selection
If a trajectory file contains valid WCBs, a `WCB_tau_*.nc` metadata file is written containing:

| Field      | Meaning                                            |
| ---------- | -------------------------------------------------- |
| `position` | Index of the trajectory inside the trajectory file |
| `ascent`   | Start time index of the ascent                     |
| `tau600`   | Duration of the ascent in time steps               |

This metadata is later used to extract full time-series variables.

The metadata is written using:
```
    write_WCB_ncfile(fl)
```
and applied over all files by:
```
    write_all_WCB_meta_files()
```
### 🔎 Why Two Cyclone Boxes?
WCBs can rise far from the cyclone center later, or be close but not ascending. The two cyclone checks ensure:
    - At `TIME_CHECK_1`: the air parcel is still in the cyclone warm sector → **identifies the correct air stream origin**
    - At `TIME_CHECK_2`: the parcel remains within the same synoptic system → **avoids false positives**

