#source("renv/activate.R")
Sys.setenv(TERM_PROGRAM = "vscode")
fpath <- file.path(Sys.getenv(
  if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"),
  ".vscode-R", "init.R")
if (file.exists(fpath)) source(fpath)
