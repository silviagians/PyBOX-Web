# ======================================================================
Box Model for polydisperse, particle laden gravity currents
# ----------------------------------------------------------------------
# This is a MODIFIED version of Pybox 1.0 adapted for web applications.
# Original authors: G. Biagioli, A. Bevilacqua, T. Esposti Ongaro, M. de' Michieli Vitturi
# Institution: Istituto Nazionale di Geofisica e Vulcanologia (INGV), Pisa, Italy
# Modified by: Silvia Giansante
# Features: box model approach, energy conoid method, automated DEM retrieval via Microsoft Planetary Computer STAC API

import numpy as np
import pandas as pd
import argparse
import sys
import scipy
import scipy.integrate
import os
import shutil
import rasterio
from pyproj import Transformer
from linecache import getline
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.transform import from_origin

# libraries for Microsoft Planetary Computer
import planetary_computer
from pystac_client import Client
import stackstac
import requests

# ========================================================================
# 1. UTILITY AND DATA LOADING FUNCTIONS
# ========================================================================
def download_dem_copernicus(lat, lon, margin_meters, output_name):
    """
    Download the Copernicus DEM (GLO-30) from Microsoft Planetary Computer.
    Data centered on given coordinates with a custom margin.
    """
    # Conversion meters to degrees
    # 1 degrees in lat ~ 111.132 meters
    margin_deg_lat = margin_meters / 111132.0 
    # Lon correction based on lat
    margin_deg_lon = margin_deg_lat / np.cos(np.radians(lat))
    bounds = (lon - margin_deg_lon, lat - margin_deg_lat, 
            lon + margin_deg_lon, lat + margin_deg_lat)

    print(f"--- Searching Copernicus DEM for Lat:{lat}, Lon:{lon} (Margin: {margin_meters})m ---")

    try:
        # Test connection with a 5s timeout
        requests.get("https://planetarycomputer.microsoft.com", timeout=5)

        catalog = Client.open(
                "https://planetarycomputer.microsoft.com/api/stac/v1",
                modifier=planetary_computer.sign_inplace,
        )
        search = catalog.search(
                collections=["cop-dem-glo-30"],
                bbox=bounds,
        )
        items = list(search.items())

        if not items:
            print("ERROR: No Copernicus DEM data found for this area.", file=sys.stderr)
            sys.exit(-11)

        # Uploads data with stackstac
        # epsg=4326 keeps Lat/Lon coordinates at the beginning
        stack = stackstac.stack(
                items, 
                assets=["data"],
                bounds_latlon=bounds, 
                epsg=4326
        )

        # Reduces data size
        dem_data = stack.mean(dim="time").sel(band="data").compute().fillna(0)

        full_path = os.path.abspath(output_name)

        # Saves DSM GeoTIFF
        with rasterio.open(
                full_path,
                'w', 
                driver='GTiff',
                height=dem_data.shape[0],
                width=dem_data.shape[1],
                count=1,
                dtype=str(dem_data.dtype),
                crs="EPSG:4326",
                transform=stack.spec.transform,
        ) as dst:
                dst.write(dem_data.values, 1)

        print(f"---Download of Copernicus DEM completed: {full_path}---")
        return full_path

    except requests.exceptions.Timeout:
        print("ERROR: Connection timeout. Data provider (Microsoft) is not responding in time.", file=sys.stderr)
        sys.exit(-13)
    except requests.exceptions.ConnectionError:
        print("ERROR: Service unavailable. Cannot reach Microsoft Planetary Computer.", file=sys.stderr)
        sys.exit(-14)
    except Exception as e:
        print(f"ERROR: Unexpected error during data retrieval: {e}", file=sys.stderr)
        sys.exit(-12)

def read_dem(dem_file, vertical_scale_factor, target_crs=None):
    """Reads DEM file (.asc and .tif) and converts it to UTM if necessary"""
    with rasterio.open(dem_file) as src:
        # Conversion from degrees to metres (UTM)
        if target_crs and src.crs != target_crs:
            print(f"Reprojecting DEM from {src.crs} to {target_crs}...")
            transform, width, height = calculate_default_transform(
                    src.crs, target_crs, src.width, src.height, *src.bounds)

            # Creates a matrix to store reprojected data
            zdem = np.zeros((height, width), dtype=np.float32)
            reproject(
                    source=rasterio.band(src, 1),
                    destination=zdem,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.bilinear)
            ncols, nrows = width, height
            cellsize = transform[0]
            xll = transform[2]
            yll = transform[5] + (nrows * transform[4]) # transform[5]=top edge, transform[4]=pixel negative dimension
        else:
            # Standard reading when data are already in UTM
            zdem = src.read(1).astype(float)
            ncols, nrows = src.width, src.height
            cellsize = src.res[0]
            xll, yll = src.bounds.left, src.bounds.bottom

        # Compute the UTM coordinates of each element of the topography from the lower-left corner
        xdem = xll + (np.arange(ncols) + 0.5) * cellsize
        ydem = yll + (np.arange(nrows) + 0.5) * cellsize
    
        # Load the DEM into a NumPy array
        zdem = np.flipud(zdem)
        zdem = np.maximum(zdem * vertical_scale_factor, 0.)

    return ncols, nrows, xll, yll, cellsize, xdem, ydem, zdem

def write_dem(output_file, ncols, nrows, xll, yll, cellsize, invasion, crs):
    """Writes the 2D invasion map as a GeoTIFF file."""
    # Creates transform (top-left origin)
    transform = from_origin(xll, yll + (nrows * cellsize), cellsize, cellsize)

    with rasterio.open(
        output_file, 'w', driver='GTiff',
        height=nrows, width=ncols,
        count=1, dtype='uint8',
        crs=crs, transform=transform,
        ) as dst:
            #Re-flip to match standard raster orientation
            dst.write(np.flipud(invasion).astype('uint8'),1)

def convert_to_utm(lat, lon):
    """Converts Geographic (Lat/Lon) coordinates to metric UTM"""
    zone = int((lon + 180) / 6) + 1 # Determine UTM zone
    proj_str = f"+proj=utm +zone={zone} +{'north' if lat >= 0 else 'south'} +datum=WGS84 +units=m +no_defs"
    transformer = Transformer.from_crs("EPSG:4326", proj_str) 
    utm_x, utm_y = transformer.transform(lat, lon)
    return utm_x, utm_y, proj_str

## Subroutine for dividing a DEM into azimuthal sectors
def dem_section(theta, nmax, xv, yv, d, ncols, nrows, cellsize, xdem, ydem, zdem, lmax):
    """Computes an altitude profile along a specific azimuth (theta)."""
    z = np.zeros(nmax + 1)
    x_pos = np.zeros(nmax + 1)
    y_pos = np.zeros(nmax + 1)
    x_pos[0], y_pos[0] = xv, yv
    theta_rad = np.radians(theta % 360)

    flag, i1 = 0, nmax
    for i in range(1, nmax + 1):
        x_pos[i] = x_pos[i-1] + d * np.cos(theta_rad)
        y_pos[i] = y_pos[i-1] + d * np.sin(theta_rad)
        if i * d > lmax:
            flag, i1 = 2, i - 1
            break
        elif (x_pos[i] > xdem[-1] + 0.5 * cellsize or
              y_pos[i] > ydem[-1] + 0.5 * cellsize or
              x_pos[i] < xdem[0] - 0.5 * cellsize or 
              y_pos[i] < ydem[0] - 0.5 * cellsize):
            flag, i1 = 1, i - 1
            break

    # Determine DEM pixel indices 

    for i in range(i1 + 1):
        ix = np.searchsorted(xdem, x_pos[i], side='left')
        iy = np.searchsorted(ydem, y_pos[i], side='left')
        # Boundary safety
        ix = np.clip(ix, 0, ncols - 1)
        iy = np.clip(iy, 0, nrows - 1)
        z[i] = zdem[iy, ix]

    return flag, z, i1

# ========================================================================
# 2. PHYSICS SUBROUTINES
# ========================================================================

def polydisperse_rg(eps, g, rg, phi):
    """Computes reduced gravity for the mixture."""
    return g * (rg + np.sum(phi * eps))

def polydisperse_density(eps, rhos, rhog):
    """Computes mixture density."""
    return np.sum(rhos * eps) + rhog * (1. - np.sum(eps))

def settling(theta, ds, rhos, rhog, kmax=100, tol=1e-4, g=9.81):
    """Calculates particle terminal velocity."""
    mu_g = 1.458e-6 * theta**(1.5) / (theta + 110.4)
    ws = ds**2 * g * rhos / (18 * mu_g)
    ws_old = ws
    Re = rhog * ds * ws / mu_g
    for _ in range(kmax):
        C_D = 24.0 / Re * (1. + 0.15 * Re**0.687)
        ws = np.sqrt(4. * rhos * g *ds / (3.* C_D *rhog))
        Re = rhog * ds * ws / mu_g
        if np.abs(ws - ws_old) < tol: break
        ws_old = ws
    return  ws if Re < 1000. else np.sqrt(4.* rhos * g * ds / (3. * 1. * rhog))

# ========================================================================
# 3. ODE FUNCTIONS FOR RKF45
# ========================================================================

def fun_cartesian(t, x, nsolid, g, rg, phi, Fr, v0, ws):
    fun_car = np.zeros(nsolid + 1)
    gprime = polydisperse_rg(x[1:], g, rg, phi)
    if gprime > 0:
        fun_car[0] = Fr * np.sqrt(gprime * v0 / x[0])
        fun_car[1:] = -ws[:] * x[1:] * x[0] / v0
    return fun_car
	
def fun_cylindrical(t, x, nsolid, g, rg, phi, Fr, v0, ws, f_geom):
    fun_cyl = np.zeros(nsolid + 1)
    gprime = polydisperse_rg(x[1:], g, rg, phi)
    if gprime > 0:
        fun_cyl[0] = Fr * np.sqrt(gprime * v0 / (f_geom * np.pi)) / x[0]
        fun_cyl[1:] = -np.pi * f_geom * ws[:] * x[1:] * (x[0]**2) / v0
    return fun_cyl
	
def event_gip(t, x, nsolid, g, rg, phi, *args):
    """Termination event: reduced gravity becomes zero."""
    return polydisperse_rg(x[1:], g, rg, phi)

# ========================================================================
# 4. MAIN BOX MODEL FUNCTION
# ========================================================================

def run_box_model(args):
    # --- Local variable setup ---
    eps0, rhos, ds = np.array(args.eps0), np.array(args.rhos), np.array(args.ds)
    nsolid = eps0.size
    
    # --- Consistency and validity checks ---

    # Check maximum number of solid phases
    Nmax = 21
    if nsolid > Nmax:
        print(f"ERROR: The number of particle classes {nsolid} exceeds the maximum limit {Nmax}.", file=sys.stderr)
        sys.exit(-1)

    # Check array dimension consistency
    if rhos.size != nsolid:
        print(f"ERROR: inconsistent number of particle densities {rhos.size} vs volume fractions {nsolid}", file=sys.stderr)
        sys.exit(-2)
    if ds.size != nsolid:
        print(f"ERROR: inconsistent number of particle sizes {ds.size} vs volume fractions {nsolid}", file=sys.stderr)
        sys.exit(-3)

    # Physical validity check: total solid volume fraction must be < 1
    total_eps = np.sum(eps0)
    if total_eps >= 1.0:
        print(f"ERROR: Total solid volume fraction ({total_eps:.4f}) is >= 1. Gas phase must be present.", file=sys.stderr)
        sys.exit(-4)

    print("Integrity checks passed. Starting simulation...\n")

    # --- Physics setup ---
    rhog = args.rhoa * 300 / args.theta0
    rg = (rhog - args.rhoa) / args.rhoa
    phi = (rhos - rhog) / args.rhoa

    # Check on initial density
    rhoc_initial = polydisperse_density(eps0, rhos, rhog)
    if rhoc_initial < args.rhoa:
        print('ERROR: Initial current density is lower than ambient density', file=sys.stderr)
        sys.exit(-5)

    # Computes settling velocities 
    ws = np.array([settling(args.theta0, d, r, rhog, g=args.g) for d, r in zip(ds, rhos)])
    # Estimates max time
    ws_min = np.min(ws)
    if int(args.theta0) != 300:
        tmax = (np.log((args.rhoa - rhog) / np.sum(eps0 * (rhos - rhog))) * (-args.h0 / ws_min))
    else:
        tmax = (np.log((np.sum(eps0) * args.rhoa) / np.sum(eps0 * rhos)))*(-args.h0 / ws_min)

    # --- DEM loading ---
    dem_data = None
    epsg_code = None
    if args.flag_DEM:
        # Calculates CRS UTM based on longitude
        zone = int((args.lon + 180) / 6) + 1
        epsg_code = f"EPSG:326{zone}" if args.lat >= 0 else f"EPSG:327{zone}"

        dem_data = read_dem(args.dem_file, args.vertical_scale_factor, target_crs=epsg_code)
        nx, ny, xll, yll, cellsize, xdem, ydem, zdem = dem_data
        if not (xdem[0] <= args.xv <= xdem[-1] and ydem[0] <= args.yv <= ydem[-1]):
            print(f"ERROR: Vent location (X:{args.xv:.2f}, Y:{args.yv}) outside DEM range.", file=sys.stderr)
            print(f"DEM Range X: {xdem[0]:.2f} to {xdem[-1]:.2f}", file=sys.stderr)
            print(f"DEM Range Y: {ydem[0]:.2f} to {ydem[-1]:.2f}", file=sys.stderr)
            sys.exit(-8)

    # --- Numerical integration ---
    y0 = np.concatenate(([args.l0], eps0))
    event_gip.terminal, event_gip.direction = True, -1

    if args.flag_coords > 0: 
        # Cartesian coordinates
        v0 = args.l0 * args.h0
        integration_fun = fun_cartesian
        extra_args = (nsolid, args.g, rg, phi, args.Fr, v0, ws)

    else:
        # Cylindrical coordinates
        spreading = (args.anglemax - args.anglemin) % 360 or 360
        f_geom = spreading / 360.
        v0 = np.pi * (args.l0**2) * args.h0 * f_geom
        integration_fun = fun_cylindrical
        extra_args = (nsolid, args.g, rg, phi, args.Fr, v0, ws, f_geom)

    # --- Numerical integration of differential equations with RKF45 method ---		
    sol = scipy.integrate.solve_ivp(
                integration_fun, [0, tmax], y0, method = 'RK45', 
                t_eval=np.arange(0, tmax, args.dt), events=event_gip, 
                args=extra_args, rtol=args.r_tol, atol=args.a_tol
                )

    t = sol.t
    l = sol.y[0] # Current front length (array)
    eps = sol.y[1:] # Concentration matrix (nsolid * len(t))
    i_max = len(t) - 1

    # Array initialisation for results
    h = np.zeros(len(t))
    u = np.zeros(len(t))
    rhom_arr = np.zeros(len(t))
    hmax = np.zeros(len(t))
    tpe = np.zeros(len(t))
    tke = np.zeros(len(t))
    deposit = np.zeros((len(t), nsolid))

    # Loop physical calculation (step by step)   
    for i in range(len(t)):
        # Calculates g' (gprime) for the polydisperse mixture
        gprime_curr = polydisperse_rg(eps[:,i], args.g, rg, phi)

        # Mixture density
        rhom_arr[i] = polydisperse_density(eps[:,i], rhos, rhog)

        # Current height based on geometry
        if args.flag_coords > 0:
            h[i] = v0 / l[i]
        else:
            h[i] = v0 / (np.pi * f_geom * (l[i]**2))

        # Front velocity
        u[i] = args.Fr * np.sqrt(np.maximum(gprime_curr * h[i], 0))

        # Energy
        tpe[i] = gprime_curr * args.rhoa * v0 * h[i]
        tke[i] = 0.5 * rhom_arr[i] * (u[i]**2)

        # Deposit matrix for _thickness.csv file
        deposit[i,:] = ws * rhos * eps[:,i] * args.dt

        # hmax: energy limit for the energy conoid
        # Avoid division per zero if rhom == rhoa
        denom = np.maximum(rhom_arr[i] - args.rhoa, 1e-6)
        hmax[i] = 0.5 * (1. / args.g) * (rhom_arr[i] / denom) * (u[i]**2)

    print("Simulation core finished. Saving files...")

    # --- Deposit & thickness calculation ---
    output_thick = args.outpfile + '_thickness.csv'

    with open(output_thick, 'w') as fID1:
        header_phases = ",".join([f"granulometric class {n} deposit thickness(m)" for n in range(nsolid)])
        fID1.write(f"current front position(m),total deposit thickness(m),{header_phases}\n")
        # Loop through time steps (from 1 to i_max)
        for i in range(len(t)):

            # Computing deposited mass from time i to final 
            # (Sum of the mass of each solid phase)
            # axis=0 sum along the time column for each solid phase
            mass_per_unit_area = np.sum(deposit[i:, :], axis=0)

            # Thickness calculation (thick = mass / (area * packing * density))
            thick = mass_per_unit_area / (args.alpha * rhos)

            # Writing front position, deposit thickness for each solid phase, total deposit thickness
            thick_str = ",".join([f"{val:3.3e}" for val in thick])
            fID1.write(f"{l[i]:5.3e},{np.sum(thick):3.3e},{thick_str}\n")

    print(f"Thickness profiles successfully saved to: {args.outpfile}_thickness.csv")

    # Writing _thickness.csv file
    with open(f"{args.outpfile}.csv", 'w') as fID2:
        header_eps = ",".join([f"eps_{n}" for n in range(nsolid)])
        fID2.write(f"length(m),height(m),rho_c(kg/m3),u(m/s),TPE(J),TKE(J),hmax(m),time(s),{header_eps}\n")
        for i in range(len(t)):
            eps_str = ",".join([f"{eps[n,i]:.6e}" for n in range(nsolid)])
            fID2.write(f"{l[i]:10.3f},{h[i]:10.3f},{rhom_arr[i]:10.3f},{u[i]:10.3e},{tpe[i]:10.3e},{tke[i]:10.3e},{hmax[i]:10.3e},{t[i]:8.2f},{eps_str}\n")
    print(f"Simulation completed. Files saved: {args.outpfile}.dat and {output_thick}")

    # --- Energy conoid (if flag_DEM is True), topography --- 
    if args.flag_DEM:
        sr = args.rad_res

        # Create the hmax(l) array interpolated for the radial distance
        l_integra = np.arange(0, int(np.max(l)) + int(sr) + 1)
        hmax_of_l = np.interp(l_integra, l, hmax)

        # Calculate invasion (2D map)
        number_of_sectors = 360
        dmax = np.zeros(number_of_sectors)

        # Define the angles to analyse
        amin = args.anglemin
        amax = args.anglemax
        if amax < amin:
            amax += 360 # To ensure anglemax > anglemin
        angles = np.arange(amin, amax + 1) % 360
        for nb in angles:
            # Call dem_section function
            _, zr, imax_sect = dem_section(nb, 2000, args.xv, args.yv, sr, nx, ny, cellsize, xdem, ydem, zdem, np.max(l))

            zrmin = zr[0]
            for na in range(1, imax_sect):
                dist_idx = int(sr * na)
                if dist_idx >= len(hmax_of_l): break

                hmx_limit = hmax_of_l[dist_idx]
                dmax[int(nb)] = sr * na

                # Energy conoid condition
                if args.differential_topography:
                    zrmin = min(zr[na], zrmin)
                    if (zr[na] - zrmin > hmx_limit): break
                else:
                    if (zr[na] > hmx_limit): break

        # Creating final 2D invasion map
        Xdem, Ydem = np.meshgrid(xdem, ydem)
        Dist = np.sqrt((Xdem - args.xv)**2 + (Ydem - args.yv)**2)
        Angle = (np.degrees(np.arctan2(Ydem - args.yv, Xdem - args.xv))) % 360

        # Boolean mask for sector and distance
        if args.anglemax < args.anglemin:
            # Exceptional case: e.g. angles from 350 to 10
            # Every angle > 350 OR < 10 is taken
            in_sector = (Angle >= args.anglemin) | (Angle <= args.anglemax)
        else:
            # Standard case: e.g. 10° to 50°
            # Each angle between 10° and 50° is considered
            in_sector = (Angle >= args.anglemin) & (Angle <= args.anglemax)
        # Mapping dmax on grid angles
        dmax_grid = dmax[Angle.astype(int)]
        # Invasion = radius shorter than limit AND inside the chosen sector
        invasion = (Dist < dmax_grid) & in_sector

        # Writing ESRII ASCII map
        write_dem(f"{args.outpfile}_EC2.tif", nx, ny, xll, yll, cellsize, invasion.astype(int), epsg_code)
        print(f"2D GeoTIFF map saved: {args.outpfile}_EC2.tif")

# ==============================================================================================
# 5. PARSER & EXECUTION FLOW
# ==============================================================================================

def parse_arguments():
    """
    Parses command-line arguments for the BoxModel simulation.
    Defines default values for all physical and numerical parameters.
    """
    parser = argparse.ArgumentParser(description="BoxModel simulation for PDCs")

    # --- 1. COMMAND LINE ARGUMENTS (to modify) ---
    parser.add_argument('--l0', type=float, default=500.0)
    parser.add_argument('--h0', type=float, default=400.0)
    parser.add_argument('--theta0', type=float, default=300.0)
    parser.add_argument('--dt', type=float, default=10.0)
    parser.add_argument('--xv', type=float, default=433015.0)
    parser.add_argument('--yv', type=float, default=4518318.99)
    parser.add_argument('--lat', type=float, default=None)
    parser.add_argument('--lon', type=float, default=None)
    parser.add_argument('--margin', type=float, default=10000) # meters

    # Parameters accepting lists (nargs='+')
    parser.add_argument('--eps0', type=float, nargs='+', default=[0.01])
    parser.add_argument('--rhos', type=float, nargs='+', default=[3000.0])
    parser.add_argument('--ds', type=float, nargs='+', default=[1e-3])
    
    # --- 2. DEFAULT CONFIGURATION PARAMETERS ---
    parser.add_argument('--Fr', type=float, default=1.2)
    parser.add_argument('--g', type=float, default=9.81)
    parser.add_argument('--rhoa', type=float, default=1.18)
    parser.add_argument('--alpha', type=float, default=0.63)
    parser.add_argument('--r_tol', type=float, default=1e-4)
    parser.add_argument('--a_tol', type=float, default=1e-6)
    parser.add_argument('--flag_coords', type=float, default=-1)
    parser.add_argument('--flag_DEM', type=float, default=1)
    parser.add_argument('--vertical_scale_factor', type=float, default=1.0)
    parser.add_argument('--rad_res', type=float, default=50.0) 
    parser.add_argument('--anglemin', type=float, default=0.0)
    parser.add_argument('--anglemax', type=float, default=360.0)
    parser.add_argument('--differential_topography', type=float, default=1)
    
    # Input/output files and cache directory
    parser.add_argument('-d', '--dem_file', type=str, default='dem.asc', help="Input DEM filename")
    parser.add_argument('-o', '--outpfile', type=str, default='boxmodeloutput', help="Output base filename")
    parser.add_argument('--cache_dir', type=str, default=None)

    return parser.parse_args()

# =======================================================================
# 6. EXECUTION FLOW
# =======================================================================
if __name__ == "__main__":
    # --- Initalises parameters ---
    args = parse_arguments()

    if args.lat is not None and args.lon is not None:
        # Converts Lat/Lon into UTM coordinates
        args.xv, args.yv, _  = convert_to_utm(args.lat, args.lon) # Easting, Northing
        print(f"Converted UTM coordinates: X={args.xv:.2f}, Y={args.yv:.2f}")

        if args.flag_DEM:
            lat_c = round(args.lat, 2)
            lon_c = round(args.lon, 2)
            margin_val = int(args.margin)

            # DEM file definition
            auto_dem_file = f"copernicus_dem_{lat_c:.2f}_{lon_c:.2f}_m{margin_val}.tif"

            # Creates 'dem_cache' subdirectory
            if args.cache_dir:
                # Uses provided path by the user
                cache_dir = os.path.abspath(args.cache_dir)
            else:
                # Default: "dem_cache" directory in the same dir as the script
                base_path = os.path.dirname(os.path.abspath(__file__))
                cache_dir = os.path.join(base_path, "dem_cache")

            # Creates the folder if absent
            if not os.path.exists(cache_dir):
                if not os.path.isdir(cache_dir):
                    os.makedirs(cache_dir)

            # Final path for cache
            full_path_cache = os.path.join(cache_dir, auto_dem_file)

            # Download DEM file if absent in the cache directory
            if not os.path.exists(full_path_cache):
                print(f"--- DEM not in cache. Downloading... ---")
                # Download data using original lat/lon
                args.dem_file = download_dem_copernicus(args.lat, args.lon, args.margin, full_path_cache)
            else:
                args.dem_file = full_path_cache
                print(f"--- DEM found in cache: {full_path_cache} ---")

            # Copy dem to working directory
            if os.path.exists(full_path_cache):
                try:
                    current_work_dir = os.getcwd()
                    new_dem_name = f"{args.outpfile}.tif"
                    destination_path = os.path.join(current_work_dir, new_dem_name)
            
                    if os.path.abspath(full_path_cache) != os.path.abspath(destination_path):
                        shutil.copy2(full_path_cache, destination_path)
                        print(f"--- DEM copied and renamed to: {destination_path} ---")

                    args.dem_file = destination_path

                except (PermissionError, OSError, shutil.Error) as e:
                    print(f"ERROR: System could not copy DEM file: {e}", file=sys.stderr)
                    sys.exit(-15) 

            else:
                print("ERROR: DEM file could not be retrieved. Skipping copy.", file=sys.stderr)
                sys.exit(-16) 

    # --- Log file creation (CSV format) ---
    log_filename = f"{args.outpfile}_params.txt"

    with open(log_filename, "w") as f:
        f.write("="*50 + "\n")
        f.write("SIMULATION INPUT LOG\n")
        f.write("="*50 + "\n")

        # Change args into a python dict for a better visualisation of input data
        for key, value in vars(args).items():
            f.write(f"{key:25} : {value}\n")
        f.write("="*50 + "\n")

    print(f"List of parameters successfully saved to: {log_filename}")

    # --- RUNNING BOXMODEL FUNCTION ---
    run_box_model(args)
