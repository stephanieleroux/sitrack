# sitrack


A Lagrangian Sea-Ice Particule Tracker for SI3 (or any sea-ice GCM running on the Arakawa C-grid).





### Dependencies

#### `mojito`
 * clone `mojito` into the `<somewhere>` directory on your computer 
 ```
 cd <somewhere>/
 git clone git@github.com:brodeau/mojito.git
 ```
 * make `python` able to locate the mojito modules, so add the following line to your `.bashrc`, `.profile`, or equivalent:
 ```
 export PYTHONPATH=<absolute_path_to_somewhere>/mojito:${PYTHONPATH}
 ```
 

### Getting started with `sitrack`

* download and extract the tarball `sitrack_demo.tar` of the directory containing the netCDF files for the demo test: `sitrack_demo`
* the absolute path to the content of this content is hereafter referred to as `<PATH_TO_DEMO_DIR>`


#### Generation of the seeding file (netCDF)

This is done by means of script `tools/generate_idealized_seeding.py` of `sitrack` 

* generate the seeding file that `si3_part_tracker.py` will use
```
generate_idealized_seeding.py -d '1996-12-15_00:00:00' \
                              -m <PATH_TO_DEMO_DIR>/NANUK4/mesh_mask_NANUK4_L31_4.2_1stLev.nc \
                              -i <PATH_TO_DEMO_DIR>/NANUK4/NANUK4-BBM23U06_1h_19961215_19970420_icemod.nc \
                              -k 0 \
                              -S 5 \
                              -f  <PATH_TO_DEMO_DIR>/NANUK4/mask_SeedInit_TrackIce_NANUK4.nc \
                              -N NANUK4
```
Basically, we are telling `generate_idealized_seeding.py`:
- the seeding occurs at `1996-12-15_00:00:00`
- find SI3 coordinates, mask, metrics, etc, into `mesh_mask_NANUK4_L31_4.2_1stLev.nc`
