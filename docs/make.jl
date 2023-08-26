using Documenter, Phaj, DocumenterTools

makedocs(sitename="P.H.A.J. Documentation")

deploydocs(
    repo = "github.com/cookienocreams/P.H.A.J."
    , branch = "gh-pages"
    , devbranch = "main"
)