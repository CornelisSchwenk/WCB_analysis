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
