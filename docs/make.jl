using Documenter, Phaj, DocumenterTools

makedocs(sitename="P.H.A.J. Documentation")

deploydocs(
    repo = "github.com/cookienocreams/Phaj"
    , branch = "docs/"
    , devbranch = "main"
)
