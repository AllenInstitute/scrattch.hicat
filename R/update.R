# Delete before releasing to public
update_iterclust <- function() {
  devtools::install_github("AllenInstitute/iterclust",
                           auth_token = "96cb6605b9e7d395b5b3e3e3c04f0eb001cf4674")
}
