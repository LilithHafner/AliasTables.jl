using OffsetTables
using Documenter

DocMeta.setdocmeta!(OffsetTables, :DocTestSetup, :(using OffsetTables); recursive=true)

makedocs(;
    modules=[OffsetTables],
    authors="Lilith Orion Hafner <lilithhafner@gmail.com> and contributors",
    sitename="OffsetTables.jl",
    format=Documenter.HTML(;
        canonical="https://LilithHafner.github.io/OffsetTables.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LilithHafner/OffsetTables.jl",
    devbranch="main",
)
