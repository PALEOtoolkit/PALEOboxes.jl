# PALEOboxes

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
CurrentModule = PALEOreactions.Forcings
```
## Forcings
```@docs
Forcings.ReactionForceInterp
GridForcings.ReactionForceGrid
FluxPerturb.ReactionFluxPerturb
FluxPerturb.ReactionRestore
```

## Grids
ReactionUnstructuredVectorGrid