# Load every BSgenome supported by this workflow before invoking ArchR. ArchR
# stores the genome as an expression in the project and evaluates it lazily in
# operations such as addGroupCoverages, so merely installing the package is not
# sufficient: its exported genome object must be attached to the R session.
supported_bsgenome_packages <- c(
  "BSgenome",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Mmusculus.UCSC.mm10",
  "BSgenome.Mmusculus.UCSC.mm39",
  "BSgenome.Rnorvegicus.UCSC.rn6"
)

missing_bsgenome_packages <- supported_bsgenome_packages[
  !vapply(
    supported_bsgenome_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]
if (length(missing_bsgenome_packages) > 0) {
  stop(
    "Missing required BSgenome package(s): ",
    paste(missing_bsgenome_packages, collapse = ", ")
  )
}

for (package_name in supported_bsgenome_packages) {
  suppressPackageStartupMessages(
    library(package_name, character.only = TRUE)
  )
}

message(
  "Loaded supported BSgenome packages: ",
  paste(supported_bsgenome_packages[-1], collapse = ", ")
)
