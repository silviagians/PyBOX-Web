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

## Original references:
* Manual: https://doi.org/10.5281/zenodo.2616551
* Scientific paper: http://dx.doi.org/10.1016/j.jvolgeores.2016.08.002

## Requirements
To run this model, you will need Python 3.x and the following libraries:
'numpy', 'pandas', 'scipy', 'rasterio', 'stackstac', 'planetary-computer'. 