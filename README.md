# PyBOX-Web: Polydisperse gravity current Box Model

This repository contains a modified version of the PyBOX 1.0 software. This adaptation is designed to function as a web-integrated application, featuring automated Digital Elevation Model (DEM) retrieval. 

## Overview
PyBOX is a box model for polydisperse, particle-laden gravity currents. This specific version extends the original code to allow:
* Automated geography: Automatic download of Copernicus DEM (GLO-30) via Microsoft Planetary Computer. 
* Coordinate conversion: Transformation from Geographic (Latitude/Longitude) to metric UTM coordinates.
* Data caching: Local storage system for DEM tiles to optimize performance and reduce API calls.

## Credits & attribution
The physical model and original logic were developed by:
Giovanni Biagioli, Andrea Bevilacqua, Tomaso Esposti Ongaro, Mattia de' Michieli Vitturi - Istituto Nazionale di Geofisica e Vulcanologia (INGV), Pisa, Italy.

Copyright (C) 2026 Silvia Giansante and PyBOX 1.0 authors

## Original references:
* Manual: https://doi.org/10.5281/zenodo.2616551
* Scientific paper: http://dx.doi.org/10.1016/j.jvolgeores.2016.08.002

## Requirements
To run this model, you will need Python 3.x and the following libraries:
'numpy', 'pandas', 'scipy', 'rasterio', 'pyproj', 'stackstac', 'pystac-client', 'requests','planetary-computer'. 

## How to use
This script is executed from the command line. You can customise the simulation by providing specific arguments.

## Key input parameters (arguments to modify)
To run a simulation, use the following syntax:
python PyBOX-Web.py --lat 40.82 --lon 14.42 --margin 5000 --outpfile vesuvius_test

| Argument            | Description                                                |
|---------------------|------------------------------------------------------------|
| --lat               | Latitude of the vent (Decimal Degrees)                     |
| --lon               | Longitude of the vent (Decimal Degrees)                    |
| --margin            | Distance from vent for DEM download (meters)               |
| --l0                | Initial front length (meters)                              |
| --h0                | Initial current height (meters)                            |
| --theta0            | Initial temperature (Kelvin)                               |
| --eps0              | Initial volume fraction of solid (list: e.g., 0.01 0.02)   |
| --rhos              | Particle density (kg/m3. List: e.g., 2000 2500)            |
| --ds                | Particle diameter (meters. List: e.g., 0.001 0.0005).      |
| --dt                | Temporal resolution of the numerical integration (seconds) |
| --outpfile          | Base name for all output files                             |

## Output description
The simulation generates five main files, using the prefix defined in --outpfile:
1. [outpfile].tif: The Digital Elevation Model (DEM) retrieved from Microsoft Planetary Computer (Copernicus GLO-30), cropped to your area of interest. 
2. [outpfile]_params.txt: A comprehensive log file containing all input parameters (physical and numerical) used for the specific run.
3. [outpfile].csv: The physical results of the simulation. It includes:
- length: Front position (m)
- height: Current thickness (m)
- rho_c: Mixture density (kg/m3)
- u: Front velocity (m/s)
- TPE/TKE: Potential and kinetic energy (J)
- hmax: Energy conoid limit (m)
- time: Simulation time (s)
4. [outpfile]_thickness.csv: Deposit thickness data, including total thickness and individual contributions for each granulometric class.
5. [outpfile]_EC2.tif: A 2D invasion map in GeoTIFF format. This represents the DEM of interest with the calculated invasion area overlaid, based on the energy conoid condition. 
