using OffsetTables
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(OffsetTables, :DocTestSetup, :(using OffsetTables); recursive=true)

makedocs(;
    authors="Lilith Orion Hafner <lilithhafner@gmail.com> and contributors",
    repo=Remotes.GitHub("LilithHafner", "OffsetTables.jl"),
    sitename="OffsetTables.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/LilithHafner/OffsetTables.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "chairmarks.lilithhafner.com"),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LilithHafner/OffsetTables.jl",
    push_preview=true,
    devbranch="main",
)
