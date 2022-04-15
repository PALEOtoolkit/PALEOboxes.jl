
# PALEOboxes.jl

PALEOboxes.jl provides the model coupler for the PALEO model framework.

PALEO provides a toolkit for constructing biogeochemical reaction-transport models of Earth system components (eg atmosphere, ocean, sediment) as a well as coupled Earth system configurations. A [YAML](https://en.wikipedia.org/wiki/YAML) format configuration file defines both the model structure (spatial domains containing biogeochemical variables and reactions that operate on them) and model parameter values, making it straightforward to add and remove biogeochemical tracers and their interactions, change the representation of eg ocean circulation, and couple components together.

The [PALEOboxes](@ref) creates the model and implements a coupler that provides a unified mechanism for
1. ‘low-level’ coupling (e.g. linking individual redox Reactions within a [`Domain`](@ref)), on which is built
2. ‘module-level’ coupling (linking e.g. atmosphere and ocean components) based on standardising names for communicating fluxes, and which enables
3. separation of biogeochemical reaction and transport. 

