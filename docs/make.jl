using BoldLab
using Documenter

DocMeta.setdocmeta!(BoldLab, :DocTestSetup, :(using BoldLab); recursive=true)

makedocs(;
    modules=[BoldLab],
    authors="jjgomezcadenas <jjgomezcadenas@gmail.com> and contributors",
    repo="https://github.com/jjgomezcadenas/BoldLab.jl/blob/{commit}{path}#{line}",
    sitename="BoldLab.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jjgomezcadenas.github.io/BoldLab.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jjgomezcadenas/BoldLab.jl",
    devbranch="main",
)
