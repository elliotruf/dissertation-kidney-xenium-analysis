# =====================================
# Project paths
# =====================================

project_root <- getwd()

project_path <- function(...) {
  file.path(project_root, ...)
}

