# Generic Reaction catalog

```@meta
CurrentModule = PALEOboxes
```

PALEOboxes includes a catalog of generic Reactions as a starting point for model construction.


## Reservoirs
```@meta
CurrentModule = PALEOboxes.Reservoirs
```
```@docs
ReactionReservoirScalar
ReactionReservoir
ReactionReservoirWellMixed
ReactionReservoirConst
ReactionReservoirForced
ReactionConst
```
### Linking to Reservoirs from a Reaction
```@docs
ReservoirLinksVector
```

## Variable Statistics
```@meta
CurrentModule = PALEOboxes.VariableStats
```
```@docs
ReactionSum
ReactionWeightedMean
ReactionAreaVolumeValInRange
```

## Fluxes
```@meta
CurrentModule = PALEOboxes
```
```@docs
Fluxes
```
```@meta
CurrentModule = PALEOboxes.Fluxes
```
```@docs
ReactionFluxTarget
ReactionFluxTransfer
```

### Contributing to flux couplers from a Reaction
```@docs
FluxContrib
```

```@meta
CurrentModule = PALEOboxes
```
## Forcings
```@docs
Forcings.ReactionForceInterp
GridForcings.ReactionForceGrid
FluxPerturb.ReactionFluxPerturb
FluxPerturb.ReactionRestore
```

## Grids
Minimal generic model grids for test purposes.

These just define the Domain size, and don't provide coordinate or metric information (cell volume, area).  Models will usually require a Reaction with a Domain-specific implementation (eg for ocean, atmosphere) that defines coordinates and volumes etc and also implements transport (advection, eddy diffusion, etc).
```@meta
CurrentModule = PALEOboxes.GridReactions
```
```@docs
ReactionUnstructuredVectorGrid
ReactionCartesianGrid
ReactionGrid2DNetCDF
```
