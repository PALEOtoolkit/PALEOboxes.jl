"""
    collate_markdown(
        io, inputpath, outputpath; 
        includes=String[], includename="doc_include.jl", imagesdirs=["images"]
    ) -> (pages::Vector, includes::Vector{String})

Recursively look through `inputpath` and subfolders for markdown files (`.md` extension),
code files `includename`, folders in `imagesdirs`, copy to `outputpath`,  return a `pages::Vector` suitable to build a tree structure
for `Documenter.jl`.    
"""
function collate_markdown(io, inputpath, outputpath; includes=String[], includename="doc_include.jl", imagesdirs=["images"])

    docs = Any[]
    filenames = readdir(inputpath)
    
    for fn in filenames
        inpath = joinpath(inputpath, fn)
        outpath = joinpath(outputpath, fn)
        if isdir(inpath)
            if fn in imagesdirs
                println(io, "cp $inpath --> $outpath")
                mkpath(outpath)
                cp(inpath, outpath; force=true)
            else
                subdocs, _ = collate_markdown(
                    io,
                    inpath,
                    outpath;
                    includes=includes,
                    imagesdirs=imagesdirs
                )
                !isempty(subdocs) && push!(docs, escape_markdown(fn)=>subdocs)
            end
        else
            if length(fn) >= 3 && fn[end-2:end] == ".md"
                push!(docs, joinpath(splitpath(outpath)[2:end])) # remove top level folder
                println(io, "cp $inpath --> $outpath")
                mkpath(outputpath)
                cp(inpath, outpath; force=true)
            elseif fn == includename
                push!(includes, inpath)
            end
        end
    end

    return (docs, includes)
end

function escape_markdown(str::AbstractString)
    str = replace(str, "_"=>"\\_")
    str = replace(str, "\$"=>"\\\$")
    str = replace(str, "*"=>"\\*")
    return str
end
