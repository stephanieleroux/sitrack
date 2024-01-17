[![DOI](https://zenodo.org/badge/622954088.svg)](https://zenodo.org/badge/latestdoi/622954088)

# sitrack

A Lagrangian Sea-Ice Particule Tracker for SI3 (or any sea-ice GCM running on the Arakawa C-grid).

### Dependencies

#### `mojito`
 * clone [`mojito`](https://github.com/brodeau/mojito) into the `<somewhere>` directory on your computer 
 ```
 cd <somewhere>/
 git clone git@github.com:brodeau/mojito.git
 ```
 * make `python` able to locate the mojito modules, so add the following line to your `.bashrc`, `.profile`, or equivalent:
 ```
 export PYTHONPATH=<absolute_path_to_somewhere>/mojito:${PYTHONPATH}
 ```

#### `python-basemap`
Required if you want to use the plotting functionality (usually triggered with `iplot=1` in the various scripts).


### Getting started with `sitrack`

* download and extract the tarball `sitrack_demo.tar` of the directory containing the netCDF files for the demo test: `sitrack_demo`
* the absolute path to the content of this directory is hereafter referred to as `<PATH_TO_DEMO_DIR>`


#### Generation of the seeding file (netCDF)

This is done by means of script `tools/generate_idealized_seeding.py` of `sitrack` 

* generate the seeding file that `si3_part_tracker.py` will use
```
generate_idealized_seeding.py -d '1996-12-15_00:00:00' \
                              -m <PATH_TO_DEMO_DIR>/NANUK4/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
                              -i <PATH_TO_DEMO_DIR>/NANUK4/NANUK4-BBM23U06_1h_19961215_19970420_icemod.nc \
                              -k 0 -S 5 \
                              -f  <PATH_TO_DEMO_DIR>/NANUK4/mask_SeedInit_TrackIce_NANUK4.nc \
                              -N NANUK4
```
Basically, we are telling `generate_idealized_seeding.py`:
- that the seeding occurs at `1996-12-15_00:00:00`
- to find SI3 coordinates, mask, metrics, etc, into `mesh_mask_NANUK4_L31_4.2_1stLev.nc`
- to use a sea-ice/open-ocean mask based on the sea-ice fraction  read in file `NANUK4-BBM23U06_1h_19961215_19970420_icemod.nc` at record `k=0`
- to use a sub-sampling of 5 in terms of model grid points
- to restrict the seeding to the horizontal subdomain defined as a mask in file `NANUK4/mask_SeedInit_TrackIce_NANUK4.nc`
- that the NEMO config is `NANUK4`

If everything goes according to plan, you should have obtained:
- an image showing you the to-be-seeded virtual buoys on a map of the Arctic: `./figs/SEEDING/sitrack_seeding_nemoTsi3_19961215_00_HSS5.png`
- the netCDF files do be used by `sitrack` that contains the location of virtual buoys to seed: `nc/sitrack_seeding_nemoTsi3_19961215_00_HSS5.nc`


Alternatively, you can create your own netCDF seeding file, provided you respect the following organization and naming convention:
```
dimensions:
        time = UNLIMITED ; // (1 currently)
        buoy = 986 ;
variables:
        int time(time) ;
                time:units = "seconds since 1970-01-01 00:00:00" ;
        int buoy(buoy) ;
        int id_buoy(buoy) ;
                id_buoy:units = "ID of buoy" ;
        float latitude(time, buoy) ;
                latitude:_FillValue = -9999.f ;
                latitude:units = "degrees north" ;
        float longitude(time, buoy) ;
                longitude:_FillValue = -9999.f ;
                longitude:units = "degrees south" ;
        float y_pos(time, buoy) ;
                y_pos:_FillValue = -9999.f ;
                y_pos:units = "km" ;
        float x_pos(time, buoy) ;
                x_pos:_FillValue = -9999.f ;
                x_pos:units = "km" ;
```

#### Tracking of the seeded virtual buoys

Now that the seeding file is generated you can fire `si3_part_tracker.py` to track the buoys.
Example:
```
si3_part_tracker.py -i <PATH_TO_DEMO_DIR>/NANUK4/NANUK4-BBM23U06_1h_19961215_19970420_icemod.nc \
                    -m <PATH_TO_DEMO_DIR>/NANUK4/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
                    -s ./nc/sitrack_seeding_nemoTsi3_19961215_00_HSS5.nc
                    -e 1997-04-20 -N NANUK4 -p 24
```
- `-s`: specifies the seeding file to use
- `-e`: specifies the end date
- `-p 24`: create an image of the positions of the buoys every 24 model records  (_i.e._ daily in this case)

Trajectories of virtual buoys are saved into file:<br>
`nc/NEMO-SI3_NANUK4_BBM23U06_tracking_nemoTsi3_idlSeed_19961215h00_19970420h00.nc`

Maps showing the positions of the buoys are generated into the `figs/tracking/` directory.




### Using non-SI3 input data: A-grid data

First, create the `coordinates_mesh_mask.nc` file based on the A-grid on which the data is provided:

```./tools/xy_arctic_to_meshmask.py -i <path_data>/20240110_hr-nersc-MODEL-nextsimf-ARC-b20240111-fv00.0.nc -o <path_data>/coordinates_mesh_mask.nc```

Keep this file.


Now, you should generate a seeding file, containing the initial position, in
space (geographical aka GPS coordinates) and time, of the buoys you wish to
track. Here, in the example we are just using the *debug* functionality that
seeds buoys with their initial positions hard-coded into function `debugSeeding()` of
`sitrack/tracking.py`:

```./tools/generate_idealized_seeding.py -d 2024-01-10_00:00:00```

The seeding file `./nc/sitrack_seeding_debug_20240110_00.nc` has been generated
and the image `./figs/SEEDING/sitrack_seeding_debug_20240110_00.png` shows you
the initial position of the buoys on the map...

Now you can run the tracking of these buoys:

    ./si3_part_tracker.py -i <path_data>/20240110_hr-nersc-MODEL-nextsimf-ARC-b20240111-fv00.0.nc \
                          -m <path_data>/coordinates_mesh_mask.nc \
                          -s ./nc/sitrack_seeding_debug_20240110_00.nc \
                          -g A -R 3 -u vxsi -v vysi -p 12
- `-g A`: input data is on the A-grid, not a C-grid as for SI3 data
- `-R 3`: nominal spatial resolution of the input data is 3 km
- `-u vxsi -v vysi`: name of `u,v` into into input netCDF file
- `-p 12`: create an image of the positions of the buoys every 12 model records (_i.e._ 12 hours in this case)

