# How to simulate in AllPix^2

## Setup
In the folder `allpixconf` is a basic configuration for simulating 1000 events with CLICPix2.

```toml
[AllPix]
log_level = "WARNING"
log_format = "DEFAULT"
number_of_events = 1000
detectors_file = "placement.conf"

[GeometryBuilderGeant4]
world_material = "air"

[DepositionPointCharge]
model = "scan"
number_of_charges = 100
source_type = "point"

[ElectricFieldReader]
model="linear"
bias_voltage=-30V
depletion_voltage=-50V

[GenericPropagation]
temperature = 293K
charge_per_step = 10

[SimpleTransfer]
max_depth_distance = 5um

[ROOTObjectWriter]
include = "DepositedCharge", "PropagatedCharge", "PixelCharge"
```

The important variables are `number_of_events` which is how many points you are going to scan. This should preferably be a number N such that N = n^3 where n is a natural number, otherwise AllPix^2 will automatically adjust to the nearest number of such type (so better be explicit then implicit). The second important variable is `number_of_charges` which is how many "test" charges you place in one volume (in allpix^2 terminology voxel) and trace them over the volume.

# Execution
To run AllPix^2 with such a configuration runs
```sh
$ allpix -c basic.conf
|17:52:02.364|  (STATUS) Welcome to Allpix^2 v1.5.0
|17:52:02.365|  (STATUS) Initialized PRNG with system entropy seed 9336901550865671610
|17:53:23.691|  (STATUS) Loaded 6 modules                     
|17:53:26.358|  (STATUS) Initialized 6 module instantiations  
|17:53:27.728|  (STATUS) Finished run of 1000 events
|17:53:27.789|  (STATUS) [F:ROOTObjectWriter] Wrote 13439 objects to 3 branches in file:
                                              /home/mpetric/AllpixDigi/allpixconf/output/data.root
|17:53:27.793|  (STATUS) Finalization completed
|17:53:27.793|  (STATUS) Executed 6 instantiations in 4 seconds, spending 65% of time in slowest instantiation GeometryBuilderGeant4
|17:53:27.794|  (STATUS) Average processing time is 4 ms/event, event generation at 244 Hz
```
