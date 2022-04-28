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

        rj = PB.create_reaction(object, PB.ReactionBase(name="test", classname="test", external_parameters=Dict{String, Any}()))
        if hasproperty(rj, :pars)
            PB.add_par(rj, rj.pars)
        end
        # println(buf, "$object  $(length(PB.get_parameters(rj))) Parameters" )
        for p in PB.get_parameters(rj)
            md = "- `$(p.name)[$(typeof(p.v))]` = "
            md *= if p.v isa AbstractString 
                "\"$(p.v)\""
            # elseif p.v isa Type || p.v isa Bool
            #     "`$(p.v)`"
            else
                "$(p.v)"
            end
            md *= isempty(p.units) ? "," : "  ($(p.units)),"
            md *= " `default_value`=$(p.default_value),"
            md *= isempty(p.allowed_values) ? "" : " `allowed_values`=$(p.allowed_values),"
            md *= " `description`=\"$(p.description)\""
            println(buf, md)
        end
        return nothing
    end


    export PARS
end
