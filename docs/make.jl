using Documenter, Phaj

makedocs(sitename="P.H.A.J. Documentation")

#makedocs(
#    sitename = "P.H.A.J. Documentation",
#    format = Documenter.HTML(),
#    modules = [Phaj]
#)

deploydocs(
    repo = "github.com/cookienocreams/P.H.A.J..git"
)