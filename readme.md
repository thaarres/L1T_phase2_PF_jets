Prerequisites: uproot, pandas, pytables, scipy, matplotlib, numpy


Import the PF candidates and AK4 jets, and run the simple 'seeded cone' reconstruction.
```
import util, algos
events, ak4jets = util.loadFromFile('debugPF.root')
jets = algos.seedConeJetsAllEvents(events)
```

Compare kinematic quantities of the newly created jet collection and the base AK4 jets:
```
import tests
compareCollections(ak4jets,jets)
```

The script *cluster.py* contains an example of loading from file, clustering a new jet collectinon and testing its performance.

Coming soon: more algos, performance analysis utilities...

### Regions
To study the impact of performing jet finding regionally, or the 'de-regionizer':
Start by doing: `import regions`
A `Region` class is provided which simply defines its edges and can check whether particles' coordinates fall within it.
One can define a `Region` and get the particles from the DataFrame which fall within it like:
```
event = events.loc[0]
region = regions.Region(0, 1, -2, -1) # a region bounded by 0 <= eta < 1 and -2 <= phi < -1
particles_in_region = event[region.in_region(event)]
```
The 36 FPGA Correlator Layer 1 regionisation of the L1 TDR is defined as a 2D array of these regions `regions.PFRegions`.

One can get a list (or list-of-lists) of particle DataFrames from a list (or list-of-lists) of regions like:
```
particles_regionized = regions.regionize(event, regions.PFRegions) # or provide a user-defined list (or list-of-lists) of regions
```

Then one can merge these regions, applying truncation, like:
```
particles_merged = regions.merge2D(particles_regionized, truncate=(32, 128), axis=0)
```
Which will limit the number of particles per group to 32 when merging the first dimension, then 128 for the final list of particles.
