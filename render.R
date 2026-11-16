if (!requireNamespace("bookdown", quietly = TRUE)) stop("Install bookdown before rendering.")
bookdown::render_book("index.Rmd", output_format = "bookdown::gitbook")
for (folder in c("scripts", "config", "metadata", "figures")) {
  if (!file.copy(folder, "_book", recursive = TRUE, overwrite = TRUE)) stop("Unable to copy resources: ", folder)
}
message("Open _book/index.html")
