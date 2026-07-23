# PyBOX-Web: Automated Box Model for particle-laden gravity current simulation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18920969.svg)](https://doi.org/10.5281/zenodo.18920969)

## Overview
PyBOX-Web is a Python implementation of the original **PyBOX** software (Biagioli et al., 2019) for simulating the propagation of gravity-driven, polydisperse pyroclastic density currents (PDCs) using a box-model approach. It predicts two-dimensional inundation areas using the energy-conoid method and topography derived from the Copernicus GLO-30 Digital Elevation Model (DEM), a Digital Surface Model (DSM)-product. 

The repository preserves the original box model formulation while introducing new functionalities aimed at simplyfing simulation setup and enabling integration into web-based applications. PyBOX-Web can be executed as a standalone command-line code or deployed as a web service. 

The main extensions include:
* **Automated topography retrieval**: automatic download of Copernicus DEM (GLO-30) data through the Microsoft Planetary Computer database. 
* **Coordinate handling**: automatic conversion from geographic (latitude/longitude) to metric UTM coordinates.
* **DSM management:** local caching of downladed tiles to improve efficiency and reduce repeated API requests.
* **Web integration**: support for browser-based workflows where users can provide model inputs and retrieve simulation outputs without directly interacting with the code.

## Credits & attribution
The original PyBOX physical and numerical models were developed by:

Giovanni Biagioli, Andrea Bevilacqua, Tomaso Esposti Ongaro, Mattia de' Michieli Vitturi.

Istituto Nazionale di Geofisica e Vulcanologia (INGV), Pisa, Italy.

* Documentation: https://doi.org/10.5281/zenodo.2616551
* Publication: http://dx.doi.org/10.1016/j.jvolgeores.2016.08.002

## Requirements
To run this model, you will need Python 3.x and the following libraries:

`numpy`, `pandas`, `scipy`, `rasterio`, `matplotlib`, `pyproj`, `stackstac`, `pystac-client`, `requests`, `planetary-computer`

## Usage
PyBOX-Web can be executed from the command line using the following syntax:

`python PyBOX-Web.py --lat 40.82 --lon 14.42 --margin 5000 --l0 500 --h0 400 --theta0 300 --eps0 0.01 --rhos 2500 --ds 0.001 --dt 5 -o vesuvius_test`

The following input parameters can be specified from the command line:

| Argument              | Description                                                                    |
|-----------------------|--------------------------------------------------------------------------------|
| `--lat`               | Vent latitude (decimal degrees)                                                |
| `--lon`               | Vent longitude (decimal degrees)                                               |
| `--margin`            | Single distance value from the vent to define the topography spatial extent (m)|
| `--l0`                | PDC front distance (m)                                                         |
| `--h0`                | PDC thickness (m)                                                              |
| `--theta0`            | PDC emperature (Kelvin)                                                        |
| `--eps0`              | Particle volume fraction                                                       |
| `--rhos`              | Particle density (kg/m³)                                                       |
| `--ds`                | Particle diameter (m)                                                          |
| `--dt`                | Numerical integration time step (s)                                            |
| `-o`                  | Base name for output files                                                     |

For simulations involving multiple particle-size classes, specify the `--eps0`, `--rhos`, and `--ds` parameters as lists, with one value per particle-size class:

`python PyBOX-Web.py --lat 40.82 --lon 14.42 --margin 5000 --l0 500 --h0 400 --theta0 300 --eps0 0.01 0.05 --rhos 2500 2700 --ds 0.001 0.0001 --dt 5 -o vesuvius_test`

Note: `--eps0`, `--rhos`, and `--ds` must contain the same number of values, with corresponding entries defining the properties of each particle-size class.

## Output
Each simulation generates the following six output files:

**1. [outpfile].tif**: input topography (GeoTIFF) based on Copernicus GLO-30 elevation data, downloaded from the Microsoft Planetary Computer. 

**2. [outpfile]_params.txt**: a text file summarising the list of physical and numerical input parameters used for the simulation.

**3. [outpfile].csv**: a table containing the simulated evolution of the PDC properties, including:

- `length`: PDC front position (m)
- `height`: PDC thickness (m)
- `rho_c`: PDC density (kg/m<sup>3</sup>)
- `u`: PDC front speed (m/s)
- `TPE/TKE`: PDC potential and kinetic energy (J)
- `hmax`: Maximum obstacle height predicted by the energy conoid approach (m)
- `time`: Simulation time step (s)

**4. [outpfile]_thickness.csv**: a table showing deposit thickness as a function of distance from the vent, including both the total deposit thickness and the contribution of each particle-size class.
  
**5. [outpfile]_EC2.tif**: PDC inundation area (GeoTIFF). The inundation extent is computed using the energy conoid approach (Orsucci, 2014; Bevilacqua, 2016) and the input topography.

**6. [outpfile].png**: a PNG image showing the simulated PDC inundation area overlaid on the input topography.

## License
GNU General Public License v3.0

Copyright (C) 2026 Silvia Giansante and PyBOX 1.0 authors


