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


#### Generation of the seeding file (netCDF)

This is done by means of script `tools/generate_idealized_seeding.py` of `sitrack` 

* generate the seeding file that `si3_part_tracker.py` will use

