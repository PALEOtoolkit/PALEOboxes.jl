"""
    DocStrings

Extends Julia `DocStringExtensions` package to provide PALEO-specific "Abbreviations" to list 
Reaction Parameters and Variables in docstrings.

exports PARS, METHODS_SETUP, METHODS_INITIALIZE, METHODS_DO

`\$(PARS)` expands to a list of Reaction Parameters
`\$(METHODS_SETUP)`, `\$(METHODS_INITIALIZE)`, `\$(METHODS_DO)` expand to a list of ReactionMethods and Variables

NB: a Reaction is created with `create_reaction`.  In addition, `\$(METHODS_DO)` etc call `register_methods(rj::AbstractReaction)`,
so may fail if this call fails with default parameters.
"""
module DocStrings
    import DocStringExtensions as DSE
    import PALEOboxes as PB

    struct Pars <: DSE.Abbreviation end

    const PARS = Pars()
    
    function DSE.format(::Pars, buf, doc)
        local docs = get(doc.data, :fields, Dict())
        local binding = doc.data[:binding]
        local object = Docs.resolve(binding)

        # println(buf, "PARS binding $binding object $object")

        try
            rj = PB.create_reaction(object, PB.ReactionBase(name="test", classname="test", external_parameters=Dict{String, Any}()))
            if hasproperty(rj, :pars)
                PB.add_par(rj, rj.pars)
            end
            # println(buf, "$object  $(length(PB.get_parameters(rj))) Parameters" )
            for p in PB.get_parameters(rj)
                md = "- `$(p.name)[$(typeof(p.v))]`="
                md *= md_value(p.v)
                md *= isempty(p.units) ? "," : "  ($(p.units)),"
                md *= " `default_value`="*md_value(p.default_value)*","
                md *= isempty(p.allowed_values) ? "" : " `allowed_values`=$(p.allowed_values),"
                md *= " `description`=\"$(escape_md(p.description))\""
                println(buf, md)
            end
        catch e
            println(buf, escape_md("PARS exception: $e"))
        end
        return nothing
    end

    struct Methods <: DSE.Abbreviation
        field::Symbol
    end

    const METHODS_SETUP = Methods(:methods_setup)
    const METHODS_INITIALIZE = Methods(:methods_initialize)
    const METHODS_DO = Methods(:methods_do)
    
    function DSE.format(methods::Methods, buf, doc)
        local docs = get(doc.data, :fields, Dict())
        local binding = doc.data[:binding]
        local object = Docs.resolve(binding)

        try
            # println(buf, "PARS binding $binding object $object")

            rj = PB.create_reaction(object, PB.ReactionBase(name="test", classname="test", external_parameters=Dict{String, Any}()))
            if hasproperty(rj, :pars)
                PB.add_par(rj, rj.pars)
            end
            
            d = PB.Domain(name="test", ID=1, parameters=Dict{String, Any}())
            rj.base.domain = d
            PB.register_methods!(rj)

            for m in getproperty(rj.base, methods.field)
                println(buf, "- `$(m.name)`")
                for v in PB.get_variables(m)
                    md = "  - "
                    md *= v.link_optional ? "\\[" : ""
                    md *= "`$(v.localname)`"
                    md *= v.link_optional ? "\\]" : ""
                    ln = PB.combine_link_name(v)
                    md *= (ln == v.localname) ? "" : " --> $(escape_md(ln))"
                    units = PB.get_attribute(v, :units)
                    md *= "  ($(units)),"
                    md *= " `$(PB.get_var_type(v))`,"
                    vfunction = PB.get_attribute(v, :vfunction)
                    md *= (ismissing(vfunction) || vfunction == PB.VF_Undefined) ? "" : " `$vfunction`,"
                    description = PB.get_attribute(v, :description)
                    md *= " `description`=\"$(escape_md(description))\""
                    
                    println(buf, md)
                end
            end
        catch e
            println(buf, escape_md("METHODS $methods exception: $e"))
        end

        return nothing
    end



    function escape_md(str::AbstractString)
        str = replace(str, "_"=>"\\_")
        str = replace(str, "\$"=>"\\\$")
        str = replace(str, "*"=>"\\*")
        return str
    end

    md_value(v) = "$v"
    md_value(v::AbstractString) = "\"$(escape_md(v))\""
    function md_value(vv::Vector{<:AbstractString})
        evv = [escape_md(v) for v in vv]
        return "$evv"
    end
    function md_value(vv::Vector)
        str = "$vv"
        str = replace(str, "["=>"\\[")
        str = replace(str, "]"=>"\\]")
        return str
    end

    export PARS, METHODS_SETUP, METHODS_INITIALIZE, METHODS_DO
end
