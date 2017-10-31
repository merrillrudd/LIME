#===============================================================
# UPDATE DESCRIPTION FILE
#===============================================================

DESCRIPTION <- readLines("DESCRIPTION")

VERSION <- as.numeric(gsub("Version: ", "", DESCRIPTION[3])) + 0.01
DESCRIPTION[3] <- paste("Version:", VERSION)

DATE <- Sys.Date()
DESCRIPTION[4] <- paste("Date:", DATE)

writeLines(DESCRIPTION, "DESCRIPTION")

# Write lsd.version()
filename <- "R/lime.version.R"
cat("#' Function to return version number\n", file = filename)
#cat("#'\n", file = filename, append = TRUE)
#cat("#' @export\n",file = filename, append = TRUE)
cat("#'\n", file = filename, append = TRUE)
cat(".onAttach <- function(libname, pkgname)\n", file = filename, append = TRUE)
cat("{\n", file = filename, append = TRUE)
cat(paste("    packageStartupMessage(\"lsd version: ", VERSION, "\n", "Compile date: ", DATE, "\n\")\n", sep = ""), file = filename, append = TRUE)
cat("}\n", file = filename, append = TRUE)

#===============================================================

roxygen2::roxygenize()
#devtools::document()
