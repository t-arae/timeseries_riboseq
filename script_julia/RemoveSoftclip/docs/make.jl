using Documenter
using RemoveSoftclip

makedocs(
    sitename = "RemoveSoftclip",
    format = Documenter.HTML(),
    modules = [RemoveSoftclip]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
