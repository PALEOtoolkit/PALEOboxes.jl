
# PALEOboxes.jl

[![CI](https://github.com/PALEOtoolkit/PALEOboxes.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PALEOtoolkit/PALEOboxes.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://PALEOtoolkit.github.io/PALEOboxes.jl/dev)

PALEOboxes.jl provides the model coupler for the PALEO model framework.

PALEO provides a toolkit for constructing biogeochemical reaction-transport models of Earth system components (eg atmosphere, ocean, sediment) as a well as coupled Earth system configurations. A [YAML](https://en.wikipedia.org/wiki/YAML) format configuration file defines both the model structure (spatial domains containing biogeochemical variables and reactions that operate on them) and model parameter values, making it straightforward to add and remove biogeochemical tracers and their interactions, change the representation of eg ocean circulation, and couple components together. Biogeochemical reaction networks can be embedded in a Fortran or C host model providing the dynamics and transport, or combined with an offline transport representation to provide a standalone model.

PALEOboxes creates the model and implements a coupler that provides a unified mechanism for
1. ‘low-level’ coupling (e.g. linking individual redox Reactions within a Domain, on which is built
2. ‘module-level’ coupling (linking e.g. atmosphere and ocean components) based on standardising names for communicating fluxes, and which enables
3. separation of biogeochemical reaction and transport. 

