suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})


#' Extract GTDB genome IDs from tree tip labels
#' Handles labels like "GB_GCA_000000000.1" or "RS_GCF_000000000.1"
#' and matches them to the `accession` column in your metadata.
#'
#' @param tree An object of class "phylo" (GTDB tree).
#' @return A tibble with columns: tip_label, genome_id (to join on).
tree_tip_ids <- function(tree) {
  tibble::tibble(genome_id = tree$tip.label)  # 1 line - USE RAW TIP.LABEL!
}
#' Check mapping between TEMPURA–GTDB accessions and GTDB tree tips
#'
#' @param tempura_tbl A table with at least `accession` column.
#' @param tree A GTDB phylo object.
#'
#' @return A list with:
#'   - n_tempura: number of rows in tempura_tbl
#'   - n_tree_tips: number of tips in tree
#'   - n_matched: number of TEMPURA rows whose accession is present as a tip
#'   - n_unmatched: number of TEMPURA rows with no matching tip
#'   - unmatched_accessions: vector of accessions not found in the tree
check_tree_mapping <- function(tempura_tbl, tree) {
  tip_ids <- tree_tip_ids(tree)

  tempura_ids <- tempura_tbl %>%
    distinct(accession) %>%
    filter(!is.na(accession))

  matched <- tempura_ids$accession[tempura_ids$accession %in% tip_ids$genome_id]
  unmatched <- tempura_ids$accession[!tempura_ids$accession %in% tip_ids$genome_id]

  res <- list(
    n_tempura         = nrow(tempura_ids),
    n_tree_tips       = length(tree$tip.label),
    n_matched         = length(matched),
    n_unmatched       = length(unmatched),
    unmatched_accessions = unmatched
  )

  cat("[TREE MAPPING CHECK]\n",
      "  TEMPURA genomes (distinct accessions): ", res$n_tempura, "\n",
      "  Tree tips:                             ", res$n_tree_tips, "\n",
      "  TEMPURA accessions matched to tips:    ", res$n_matched, "\n",
      "  TEMPURA accessions NOT in tree:        ", res$n_unmatched, "\n",
      sep = "")

  res
}

#' Prune GTDB tree to TEMPURA–GTDB genomes
#'
#' @param tempura_tbl A table with an `accession` column (e.g. archaea_key).
#' @param tree A GTDB phylo object.
#'
#' @return A pruned phylo tree containing only TEMPURA accessions
#'         that are present as tips.
prune_tree_to_tempura <- function(tempura_tbl, tree) {
  mapping_check <- check_tree_mapping(tempura_tbl, tree)
  
  if (mapping_check$n_matched == 0) {
    stop("No TEMPURA accessions found in tree. Check accession formats.")
  }
  
  keep_ids <- tempura_tbl %>%
    distinct(accession) %>%
    filter(!is.na(accession), accession %in% tree$tip.label)
  
  pruned <- ape::keep.tip(tree, keep_ids$accession)
  
  cat("[TREE PRUNING]\n",
      "  Original tips:", length(tree$tip.label), "\n",
      "  TEMPURA matched:", mapping_check$n_matched, "\n", 
      "  Pruned tips:", length(pruned$tip.label), "\n",
      sep = "")
  
  pruned
}

#' Analyze GTDB-TEMPURA matches vs tree coverage
#' 
#' @param gtdb_tempura_file Path to GTDB-TEMPURA match file (ar_ALL_MATCH.csv or bac_ALL_MATCH.csv)
#' @param tempura_key Tempura key table (archaea_key or bacteria_key) 
#' @param tree phylo object (ar_tree or bac_tree)
#' @param tree_name String for reporting ("Archaea ar53" or "Bacteria bac120")
#' 
#' @return List with mapping results, summary table, and details
report_tree_coverage <- function(gtdb_tempura_file, tempura_key, tree, tree_name) {
  
  gtdb_tempura_ar <- read_csv(gtdb_tempura_file, show_col_types = FALSE)
  cat("\n[", tree_name, "GTDB-TEMPURA Analysis]\n")
  cat("Loaded:", nrow(gtdb_tempura_ar), "rows,", n_distinct(gtdb_tempura_ar$accession), "unique GTDB accessions\n")
  
  # 1. TOTAL GTDB vs tree
  gtdb_accessions <- unique(gtdb_tempura_ar$accession)
  gtdb_key <- tibble(accession = gtdb_accessions)
  total_mapping <- check_tree_mapping(gtdb_key, tree)
  
  # 2. Direct TEMPURA vs tree  
  direct_mapping <- check_tree_mapping(tempura_key, tree)
  missing_ar <- direct_mapping$unmatched_accessions
  
  # 3. Missing species → GTDB alts vs tree
  missing_genera <- tempura_key %>%
    filter(accession %in% missing_ar) %>%
    pull(genus_and_species) %>%
    unique()
  
  cat(length(missing_ar), "tree-missing TEMPURA →", length(missing_genera), "unique genus_species\n")
  # 4. Find ALL GTDB accessions for those genus_species
  genus_gtdb_matches <- gtdb_tempura_ar %>%
    filter(genus_and_species %in% missing_genera) %>%  
    count(genus_and_species, accession, sort = TRUE) %>%
    mutate(in_tree = accession %in% tree$tip.label)
    
  genus_gtdb_accessions <- unique(genus_gtdb_matches$accession)
  gtdb_key <- tibble(accession = genus_gtdb_accessions)
  missing_mapping <- check_tree_mapping(gtdb_key, tree)
  
  # PERFECT SUMMARY (no tibble errors!)
  cat("\n🎉 SUMMARY", tree_name, ":\n")
  cat("• Total GTDB accessions:", n_distinct(gtdb_tempura_ar$accession),
    "→", total_mapping$n_matched, 
    round(total_mapping$n_matched/n_distinct(gtdb_tempura_ar$accession)*100, 1),
    "% in tree\n")
  cat("• Direct TEMPURA-tree:", direct_mapping$n_matched, "/", nrow(tempura_key),
    "=", round(direct_mapping$n_matched/nrow(tempura_key)*100, 0), "%\n") 
  cat("•", length(missing_ar), "missing →", length(genus_gtdb_accessions), 
    "GTDB alts →", missing_mapping$n_matched, "=", 
    round(missing_mapping$n_matched/length(genus_gtdb_accessions)*100, 1),
    "% tree coverage\n")  
  list(
    total = total_mapping,
    direct = direct_mapping, 
    missing = missing_mapping,
    details = list(missing_genera = missing_genera, genus_matches = genus_gtdb_matches)
  )
}