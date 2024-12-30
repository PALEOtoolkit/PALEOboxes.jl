module CommonDataModelExt

import PALEOboxes as PB
import CommonDataModel as CDM

#####################################
# Wrapper types
######################################

struct ModelCDM <: CDM.AbstractDataset
    model::PB.Model
    modeldata::Union{Nothing, PB.AbstractModelData}
end

PB.CDModel(model::PB.Model) = ModelCDM(model, nothing)
PB.CDModel(modeldata::PB.AbstractModelData) = ModelCDM(modeldata.model, modeldata)

struct DomainCDM <: CDM.AbstractDataset
    domain::PB.Domain
    modeldata::Union{Nothing, PB.AbstractModelData}
end

PB.CDModel(domain::PB.Domain) = DomainCDM(domain, nothing)


struct VariableDomainCDM <: CDM.AbstractDataset
    variabledomain::PB.VariableDomain
    modeldata::Union{Nothing, PB.AbstractModelData}
    data
end

PB.CDModel(variabledomain::PB.VariableDomain) = VariableDomainCDM(variabledomain, nothing, nothing)

######################################
# Model
######################################

# iterable with all group names
CDM.groupnames(m::ModelCDM) = [d.name for d in m.model.domains]

CDM.group(m::ModelCDM, name::AbstractString) = DomainCDM(PB.get_domain(m.model, name; allow_not_found=false), m.modeldata)

###############################################
# Domain
################################################

CDM.name(d::DomainCDM) = d.domain.name

# TODO
# parentdataset(d::DomainCDM)

# returns a list of variable names as strings
Base.keys(d::DomainCDM) = [v.name for v in PB.get_variables(d.domain)]

function CDM.variable(d::DomainCDM, varname::AbstractString)
    variabledomain = PB.get_variable(d.domain, varname; allow_not_found=false)
    if isnothing(d.modeldata)
        data = nothing
    else
        data = PB.get_data(variabledomain, d.modeldata)
    end
    return VariableDomainCDM(variabledomain, d.modeldata, data)
end

CDM.dimnames(d::DomainCDM) = [nd.name for nd in PB.get_dimensions(d.domain)]

CDM.dim(d::DomainCDM, name::AbstractString) = PB.get_dimension(d.domain, name).size


###############################################
# VariableDomain
################################################

CDM.name(v::VariableDomainCDM) = v.variabledomain.name

CDM.dataset(v::VariableDomainCDM) = DomainCDM(v.variabledomain.domain, v.modeldata)

CDM.dimnames(v::VariableDomainCDM) = [nd.name for nd in PB.get_dimensions(v.variabledomain)]

Base.ndims(v::VariableDomainCDM) = length(CDM.dimnames(v))

Base.size(v::VariableDomainCDM) = (nd.size for nd in PB.get_dimensions(v.variabledomain))

CDM.attribnames(v::VariableDomainCDM) = keys(v.variabledomain.attributes)

CDM.attrib(v::VariableDomainCDM, name::Symbol) = v.variabledomain.attributes[name]

Base.getindex(v::VariableDomainCDM, indices...) = Base.getindex(v.data, indices...)

Base.eltype(v::VariableDomainCDM) = Base.eltype(v.data)


end # module