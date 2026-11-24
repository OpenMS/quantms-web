# Which protein modifications (PTMs, as UniProt annotates them) are enriched at
# the top or the bottom of a ranked protein list: protein set enrichment
# analysis (PSEA) with the PEIMAN2 R package.
#
# Rscript psea.R <output_dir> <input_csv> <organism> <permutations>
#   input_csv: two columns, the UniProt accession and the protein's score
#              (the log2 fold change of the group comparison)
#   organism:  as UniProt names it, e.g. "Homo sapiens (Human)"
# Writes <output_dir>/psea_results.csv (one row per modification) and
# psea_results.rds (everything runPSEA returns).
#
# Written by Payman Nickchi, the author of PEIMAN2 (the organism name is
# resolved against PEIMAN2's table before the run, see resolve_organism).

suppressPackageStartupMessages(library(PEIMAN2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: Rscript psea.R <output_dir> <input_csv> <organism> <permutations>")
}
output_dir <- args[[1]]
input_file <- args[[2]]
organism <- args[[3]]
nperm <- suppressWarnings(as.integer(args[[4]]))

if (is.na(nperm) || nperm < 1) {
    stop("The number of permutations must be a positive integer.")
}
if (!nzchar(organism)) {
    stop("The organism cannot be empty.")
}

protein <- read.csv(input_file, stringsAsFactors = FALSE)
if (ncol(protein) != 2) {
    stop("The input must have two columns: the UniProt accession and the score.")
}
protein[[2]] <- suppressWarnings(as.numeric(protein[[2]]))
if (nrow(protein) == 0 || anyNA(protein) || any(!is.finite(protein[[2]])) ||
    any(trimws(protein[[1]]) == "")) {
    stop("The input contains missing, invalid or empty values.")
}

# PEIMAN2 keeps the proteins whose organism is exactly os.name (the OS column of
# its table): resolve the name given to one the table holds. A name that differs
# only in case, full stops or the parts in brackets (UniProt FASTA headers write
# "Homo sapiens" for "Homo sapiens (Human)") takes the matching one; a name that
# matches none takes the organism most of the input proteins belong to.
resolve_organism <- function(organism, accessions) {
    known <- PEIMAN2:::load_peiman_database("bundled")
    names <- unique(known$OS[!is.na(known$OS)])
    if (organism %in% names) {
        return(organism)
    }
    scientific <- function(x) tolower(trimws(gsub(".", "", sub("\\s*\\(.*$", "", x), fixed = TRUE)))
    of_input <- table(known$OS[known$AC %in% accessions])
    same <- names[which(scientific(names) == scientific(organism))]
    if (length(same) > 0) {
        found <- of_input[names(of_input) %in% same]
        pick <- if (length(found) > 0) names(found)[which.max(found)] else same[[1]]
        cat("Organism", shQuote(organism), "is", shQuote(pick), "in PEIMAN2\n")
        return(pick)
    }
    if (length(of_input) > 0) {
        pick <- names(of_input)[which.max(of_input)]
        cat("PEIMAN2 does not know", shQuote(organism), "- most proteins belong to", shQuote(pick), "\n")
        return(pick)
    }
    stop("PEIMAN2 knows neither the organism ", shQuote(organism), " nor the proteins.")
}
organism <- resolve_organism(organism, trimws(protein[[1]]))

cat("Organism:", organism, "\n")
cat("Permutations:", nperm, "\n")
cat("Proteins:", nrow(protein), "\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
psea_res <- runPSEA(protein = protein, os.name = organism, nperm = nperm)

write.csv(psea_res[[1]], file.path(output_dir, "psea_results.csv"), row.names = FALSE)
saveRDS(psea_res, file.path(output_dir, "psea_results.rds"))
cat("PTM enrichment finished\n")
