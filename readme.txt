# `Burgers`: Numerical Schemes for Advection-Diffusion Equations

This repository contains the source code for the Mini-Project 1 for MAE557: Simulation and Modeling of Fluid Flows. The project involves implementing, comparing, and analyzing various numerical schemes for solving the one-dimensional linear and non-linear (Burgers') advection-diffusion equations.

The implemented schemes include:
*   **Temporal Discretization:** Explicit (Forward) Euler and Implicit (Backward) Euler.
*   **Spatial Discretization:** First-Order Upwind and Second-Order Centered Difference.

## Requirements

To compile and run this project, you will need the following dependencies installed:

*   **C++ Compiler:** A modern C++ compiler (the script is configured for `clang++`, but `g++` can be substituted).
*   **Armadillo:** A C++ library for linear algebra.
    *   **Note:** The `run.sh` script is configured for a Homebrew installation on macOS (`/opt/homebrew/`). Users on other systems (e.g., Linux) will need to update the library and include paths.
*   **Python 3:** Required for running the analysis and plotting scripts.
*   **Python Libraries:** `numpy` and `matplotlib`.
    *   These can be installed via pip: `pip install numpy matplotlib`
*   **FFmpeg:** Required for creating an MP4 video from the output plots.
    *   **Note:** The script uses a hardcoded path to FFmpeg. This will likely need to be adjusted.

## How to Run a Simulation

The primary method for running a simulation is through the `run.sh` script, which handles compilation, execution, and output generation.

### 1. Configure the Simulation in `run.sh`

Open the `run.sh` script in a text editor. The simulation parameters are defined in the variables at the top of the file.

| Variable | Description | Example Value |
| :--- | :--- | :--- |
| `RUNNAME` | A descriptive name for the simulation. This will also be the name of the output directory. | `"implicit_upwind_n128"` |
| `TIMEINT` | Selects the time integration scheme. **0** for Explicit Euler, **1** for Implicit Euler. | `1` |
| `SPATINT` | Selects the spatial discretization scheme. **0** for Centered Difference, **1** for Upwind. | `1` |
| `DIFFCOEF`| The diffusion coefficient ($\nu$) for the equation. | `0.01` |
| `NX1` | The number of grid points for the spatial domain. | `128` |
| `MOVIEFPS`| The frames per second for the output MP4 video. | `20` |

### 2. Platform-Specific Configuration

The `run.sh` script contains hardcoded paths that may need to be changed for your system.

*   **Compiler and Armadillo Paths:** The `clang++` command specifies paths for the Armadillo library installed via Homebrew on Apple Silicon.
    ```bash
    # Original command in script
    clang++ -I/opt/homebrew/include -L/opt/homebrew/lib -larmadillo ...

    # Example for a standard Linux installation
    # g++ -std=c++11 ../src/main.cpp -o burger -larmadillo
    ```
*   **FFmpeg Path:** The path to the FFmpeg executable is also hardcoded.
    ```bash
    # Original command in script
    /opt/homebrew/Cellar/ffmpeg/5.1.2/bin/ffmpeg -r ${MOVIEFPS} ...

    # On most systems, you can simply use the command directly:
    # ffmpeg -r ${MOVIEFPS} ...
    ```

### 3. Execute the Script

Once configured, run the script from the project's root directory:
```bash
bash run.sh
```

The script will create a directory for the run, compile the code, execute the simulation, generate plots, create a video, and organize all output files.

## Generating Analysis Plots

The repository includes separate scripts for generating the stability and convergence plots featured in the project report.

### Stability Plots

To generate the stability diagrams:
1.  Navigate to the `stability/` directory.
2.  Open the `stabilityplots.py` script in a text editor to set the desired values for the convection (`c`) and diffusion (`nu`) coefficients.
3.  Run the script from the command line:
    ```bash
    python stabilityplots.py
    ```

### Convergence Plots

To generate the convergence plot (`L2` norm of the error vs. grid spacing):
1.  Navigate to the `convergence/` directory.
2.  The script `convergenceplots.py` reads pre-generated simulation data from the `conv_data/` subdirectory.
3.  Run the script:
    ```bash
    python convergenceplots.py
    ```

## Code Structure

The project is organized into the following directories and files:

```
.
├── run.sh                # Main script to configure, compile, and run a simulation
├── report.pdf            # Project report
│
├── src/                  # Contains all C++ source and Python plotting scripts
│   ├── main.cpp          # Main C++ driver, parses arguments and runs the simulation
│   ├── grid.hpp          # Header for grid generation and management
│   ├── implicit.hpp      # Header for the implicit solver implementation
│   └── plotter.py        # Python script to generate plots from simulation output
│
├── stability/            # Scripts for generating stability plots
│   └── stabilityplots.py
│
├── convergence/          # Scripts and data for generating convergence plots
│   ├── convergenceplots.py
│   └── conv_data/
│       └── ... (CSV files with pre-computed error data)
│
└── movies/               # Contains animations used for debugging during early development
```
