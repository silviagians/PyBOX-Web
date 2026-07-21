# PyBOX-Web: Automated Box Model for particle-laden gravity current simulation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18920969.svg)](https://doi.org/10.5281/zenodo.18920969)

## Overview
PyBOX-Web is a Python implementation of the original **PyBOX** software (Biagioli et al., 2019) for the simulation of polydisperse pyroclastic density currents (PDCs). 

The repository preserves the original box model formulation while introducing new functionalities aimed at simplyfing simulation setup and enabling integration into web-based applications. PyBOX-Web can be executed as a standalone command-line code or deployed as a web service. 

The main extensions include:
* **Automated topography retrieval**: automatic download of Copernicus DEM (GLO-30) data through the Microsoft Planetary Computer database. 
* **Coordinate handling**: automatic conversion from geographic (latitude/longitude) to metric UTM coordinates.
* **DSM management:** local caching of downladed tiles to improve efficiency and reduce repeated API requests.
* **Web integration**: support for browser-based workflows where users can provide model inputs and retrieve simulation outputs without directly interacting with the code.

## Credits & attribution
The original PyBOX physical and numerical models were developed by:

Giovanni Biagioli, Andrea Bevilacqua, Tomaso Esposti Ongaro, Mattia de' Michieli Vitturi
Istituto Nazionale di Geofisica e Vulcanologia (INGV), Pisa, Italy.

PyBOX-Web builds upon this original code by extending the software architecture and adding new functionalities for automated data acquisition, preprocessing, and execution.

## Original references:
* Documentation: https://doi.org/10.5281/zenodo.2616551
* Publication: http://dx.doi.org/10.1016/j.jvolgeores.2016.08.002

## Requirements
To run this model, you will need Python 3.x and the following libraries:

`numpy`, `pandas`, `scipy`, `rasterio`, `pyproj`, `stackstac`, `pystac-client`, `requests`, `planetary-computer`

## Usage
PyBOX-Web can be executed from the command line using the following syntax:

`python PyBOX-Web.py --lat 40.82 --lon 14.42 --margin 5000 --l0 500 --h0 400 --theta0 300 --eps0 0.01 --rhos 2500 --ds 0.001 --dt 5 -o vesuvius_test`

For simulations involving multiple particle size classes, provide the corresponding parameters as lists:

`python PyBOX-Web.py --lat 40.82 --lon 14.42 --margin 5000 --l0 500 --h0 400 --theta0 300 --eps0 0.01 0.05 --rhos 2500 2700 --ds 0.001 0.0001 --dt 5 -o vesuvius_test`

The following input parameters can be specified from the command line:

| Argument              | Description                                                |
|-----------------------|------------------------------------------------------------|
| `--lat`               | Latitude of the vent (decimal Degrees)                     |
| `--lon`               | Longitude of the vent (decimal Degrees)                    |
| `--margin`            | Distance from vent used for DSM download (m)               |
| `--l0`                | Initial PDC front length (m)                               |
| `--h0`                | Initial PDC thickness (m)                                  |
| `--theta0`            | Initial temperature (Kelvin)                               |
| `--eps0`              | Initial solid volume fraction                              |
| `--rhos`              | Particle density (kg/m³).                                 |
| `--ds`                | Particle size class (m)                                    |
| `--dt`                | Numerical integration time step (s)                        |
| `-o`                  | Base name for output files                                 |

Note: For polydisperse simulations, `--eps0`, `--rhos`, and `--ds` must contain the same number of values, corresponding to the different particle-size classes.

## Output
Each simulation generates the following output files:

**1. [outpfile].tif**: Digital Surface Model (DSM) downloaded from Microsoft Planetary Computer (Copernicus GLO-30) and cropped to the selected area of interest. 

**2. [outpfile]_params.txt**: A log file containing all physical and numerical input parameters used for the simulation run.

**3. [outpfile].csv**: Time-dependent simulation outputs, including:

- `length`: PDC front position (m)
- `height`: PDC thickness (m)
- `rho_c`: PDC density (kg/m<sup>3</sup>)
- `u`: PDC front speed (m/s)
- `TPE/TKE`: PDC potential and kinetic energy (J)
- `hmax`: Maximum obstacle height predicted by the energy conoid approach (m)
- `time`: Simulation time (s)

**4. [outpfile]_thickness.csv**: Deposit thickness data, including total deposit thickness and the contribution of each particle class.
  
**5. [outpfile]_EC2.tif**: GeoTIFF map showing the simulated invasion area. The calculated inundation extent is overlaid on the input DSM according to the energy conoid approach (Orsucci, 2014; Bevilacqua, 2016).

## License
GNU General Public License v3.0

Copyright (C) 2026 Silvia Giansante and PyBOX 1.0 authors


