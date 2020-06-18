# AllPix^2 Digitiser

## Simulation preparation

Before using the digitiser the simulation has to be setup in the proper way. This refers to the stepping in sensitive silicon detector.
To achieve this three elements need to be configured

1. Step length in the region: For this one needs to define in the compact file a `limitset`
  ```XML
  <limits>
    <limitset name="allpix_digi_limit">
      <limit name="step_length_max" particles="*" value="1" unit="micrometer" />
    </limitset>
  </limits>
  ```
  Here a maximum step of 1 um is defined, this is configurable parameter and depends on how you have configured you AllPix^2 simulation and the output you have extracted (mesh size).
  You have to append this `limitset` to a detector, in the case of `ZPlanarTracker` the limits have to be attached to the `detector` tag.
  ```XML  
  <detector name="VertexBarrel" type="ZPlanarTracker" vis="VXDVis" id="DetID_VXD_Barrel" readout="VertexBarrelCollection"  region="VertexBarrelRegion" limits="allpix_digi_limit">
  ```
  The exact position of the `limits` attribute depends on the implementation of the individual detector.
  Some detector have limits for the whole detector some have for individual volumes (a preferable implementation is latter as you only want to have small stepping in the sensitive volume).

2. Sensitive action needs to store each individual step: For this one has to set the tracker action to the following
  ```py
  SIM.action.tracker=('Geant4TrackerWeightedAction', {'HitPositionCombination': 2, 'CollectSingleDeposits': True})
  ```
  where `HitPositionCombination: 2` will set hit position between the pre and post step endpoint and `CollectSingleDeposits: True` will write out each step.
  Other tracker custom tracker actions can be used if the same above mentioned boundary conditions are fulfilled.

3. Configure filter to retain all particles: As DD4hep/DDSim applies filters to simulated particles it is important to disable this filter for energy deposition. A simple configuration with 0 energy cuts is
  ```py
  ##  list of filter objects: map between name and parameter dictionary
  SIM.filter.filters = {'edep0': {'parameter': {'Cut': 0.0}, 'name': 'EnergyDepositMinimumCut/Cut0'}, 'geantino': {'parameter': {}, 'name': 'GeantinoRejectFilter/GeantinoRejector'}}
  ##  default filter for tracking sensitive detectors; this is applied if no other filter is used for a tracker
  SIM.filter.tracker = "edep0"
  ##  default filter for calorimeter sensitive detectors; this is applied if no other filter is used for a calorimeter
  SIM.filter.calo = "edep0"
  ```
  Beware there is a default 1 keV energy cut in `ddsim` if you do not explicitly disable this as shown above.

**Caveat**: To create a electron-hole pair you need 3.62 eV, all the steps that are below this value are disregarded in the digitisation. If you setup a too low stepping there might be a disproportionate amount of such steps. You might wish to directly apply a filter in DD4hep/ddsim for 3.62 eV to minimise the amount of storage.
