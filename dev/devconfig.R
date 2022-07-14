## add dependecies
# pak::pkg_install("attachment")
attachment::att_amend_desc()

## update document
devtools::document(".")

## load all
devtools::load_all(".")

## add readme
# usethis::use_readme_rmd()

## update readme
devtools::build_readme()

## more document
# usethis::use_pkgdown()
# usethis::use_pkgdown_github_pages()
# use_github_action("pkgdown")
