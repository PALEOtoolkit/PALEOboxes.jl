


#################################################
# Dimensions
#####################################################

"""
    NamedDimension(name, size)

A named dimension
"""
struct NamedDimension
    name::String
    size::Int64
end

function Base.show(io::IO, nd::NamedDimension)
    print(io, "NamedDimension(name=", nd.name, ", size=", nd.size, ")")
    return nothing
end

"""
    function get_dimensions(obj) -> Vector{NamedDimension}

Get all dimensions for PALEO object `obj`
"""
function get_dimensions end

"""
    function get_dimension(obj, dimname) -> NamedDimension

Get all dimension `dimname` for PALEO object `obj`
"""
function get_dimension end

"""
    function set_coordinates!(obj, dimname, coordinates::Vector{String})

Set coordinates (variable names) attached to `dimname` for PALEO object `obj`

PALEO convention is that where possible `coordinates` contains:
- three variable names, for cell midpoints, lower edges, upper edges, in that order.
- two variable names, for cell midpoints and a bounds variable (with a bounds dimension as the first dimensions)
"""
function set_coordinates! end

"""
    function get_coordinates(obj, dimname) -> coordinates::Vector{String}

Get coordinates (if any) attached to `dimname` for PALEO object `obj`

PALEO convention is that where possible `coordinates` contains:
- three variable names, for cell midpoints, lower edges, upper edges, in that order.
- two variable names, for cell midpoints and a bounds variable (with a bounds dimension as the first dimensions)
"""
function get_coordinates end
