# analysis/functions/metadata_integration.R
# Map TEMPURA species to GTDB genomes (Archaea and Bacteria) using
# NCBI TaxID (primary) and genus+species names (secondary).

# TEMPURA columns used:
#   row_id, genus_and_species, taxonomy_id, strain, superkingdom, ...
# GTDB metadata columns used:
#   accession, gtdb_taxonomy, ncbi_taxid, ncbi_organism_name,
#   genome_size, gc_percentage, ...

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# Small helper to normalize species names
normalize_species_name <- function(x) {
  x %>%
    str_trim() %>%
    str_squish() %>%
    str_replace_all("_", " ") %>%
    str_to_lower()
}
# Take first two words (genus + species) from a normalized name
first_two_words <- function(x) {
  stringr::word(x, 1, 2) 
}

#' Map a TEMPURA table (one superkingdom) to GTDB metadata
#'   using NCBI TaxID first, then genus+species names.
#'
#' @param tempura_tbl TEMPURA subset (Archaea or Bacteria) with:
#'   - row_id
#'   - taxonomy_id
#'   - genus_and_species
#'   - superkingdom
#' @param gtdb_meta_tbl GTDB metadata table (ar53_metadata_r226.tsv or bac120_metadata_r226.tsv) with:
#'   - accession
#'   - ncbi_taxid
#'   - ncbi_organism_name
#'   - gtdb_taxonomy
#'   - gtdb_representative
#'
#' @return List with:
#'   - mapped: TEMPURA rows joined to GTDB genomes via TaxID or name
#'   - unmapped: TEMPURA rows with no GTDB match by either route
#'   - summary: list of counts
#'
map_tempura_to_gtdb <- function(tempura_tbl, gtdb_meta_tbl) {
  required_cols_tmp  <- c("row_id", "taxonomy_id", "genus_and_species", "superkingdom")
  missing_tmp <- setdiff(required_cols_tmp, names(tempura_tbl))
  if (length(missing_tmp) > 0) {
    stop("tempura_tbl is missing required columns: ", paste(missing_tmp, collapse = ", "))
  }

  required_cols_gtdb <- c("accession", "ncbi_taxid", "ncbi_organism_name",
                         "gtdb_taxonomy", "gtdb_representative")
  missing_gtdb <- setdiff(required_cols_gtdb, names(gtdb_meta_tbl))
  if (length(missing_gtdb) > 0) {
    stop("gtdb_meta_tbl is missing required columns: ", paste(missing_gtdb, collapse = ", "))
  }

  # ------------------------------------------------------------------
  # 1. Primary join: NCBI TaxID
  # ------------------------------------------------------------------
  joined_taxid <- tempura_tbl %>%
    left_join(
        gtdb_meta_tbl, 
        by = c("taxonomy_id" = "ncbi_taxid"),
        relationship = "many-to-many", 
        multiple = "all") %>%
    select(all_of(c(names(tempura_tbl))), everything()) %>%
    mutate(
        taxid_match = !is.na(accession),
        name_only_match = FALSE
    ) %>%
    mutate(  # SECOND mutate - NOW gtdb_representative exists
        gtdb_representative = {
            cat( "taxid ncols", paste(length(names(.))), "\n")
            coalesce(gtdb_representative, FALSE)
        }
    )

  # TEMPURA rows with any TaxID match
  mapped_taxid <- joined_taxid %>%
    filter(taxid_match)

  # TEMPURA rows that are still unmapped after TaxID join
  still_unmapped_tempura <- joined_taxid %>%
    filter(!taxid_match) %>%
    distinct(row_id, .keep_all = TRUE)

  # ------------------------------------------------------------------
  # 2. Secondary join: genus_and_species ↔ first two words of
  #    ncbi_organism_name (binomial), only for rows still unmapped
  #    by TaxID.
  # ------------------------------------------------------------------
  if (nrow(still_unmapped_tempura) > 0) {
    tmp_names <- still_unmapped_tempura %>%
      select(all_of(names(tempura_tbl))) %>% 
      mutate(
        tempura_name_norm = normalize_species_name(genus_and_species),
        tempura_binomial  = first_two_words(tempura_name_norm)
        )

    gtdb_names <- gtdb_meta_tbl %>%
      #select(all_of(required_cols_gtdb)) %>%  # Includes gtdb_representative
      mutate(
        gtdb_name_norm = normalize_species_name(ncbi_organism_name),
        gtdb_binomial  = first_two_words(gtdb_name_norm)
      )

    joined_name <- tmp_names %>%
      left_join(
        gtdb_names,
        by = c("tempura_binomial" = "gtdb_binomial"),
        relationship = "many-to-many",
        multiple = "all"
        ) %>%
      select(all_of(c(names(tempura_tbl))), everything()) %>%
      mutate(
        taxid_match     = FALSE,
        name_only_match = !is.na(accession)
        ) %>%
      mutate(
        gtdb_representative = {
          cat( "names ncols", paste(length(names(.))), "\n")
          coalesce(gtdb_representative, FALSE)
      }
    )

    mapped_name <- joined_name %>%
      filter(name_only_match)

  } else {
    mapped_name <- tempura_tbl[0, ]
  }
  # ------------------------------------------------------------------
  # 3. Combine mapped sets and compute unmapped rows
  # ------------------------------------------------------------------
  mapped_combined <- bind_rows(
    mapped_taxid,
    mapped_name
  )

  # VERIFICATION CHECKS (REMOVE AFTER CONFIRMING 0 NAs)
  #cat("Topt_ave NAs:", sum(is.na(mapped_combined$Topt_ave)), "\n")
  #cat("accession NAs:", sum(is.na(mapped_combined$accession)), "\n")
  #cat("gtdb_representative NAs:", sum(is.na(mapped_combined$gtdb_representative)), "\n")

  # Which TEMPURA row_ids have any match?
  matched_row_ids <- mapped_combined %>%
    distinct(row_id) %>%
    pull(row_id)

  unmapped_final <- tempura_tbl %>%
    filter(!row_id %in% matched_row_ids)

  # Counts
  n_tempura_total <- nrow(tempura_tbl)
  n_tempura_taxid_mapped <- mapped_taxid %>% distinct(row_id) %>% nrow()
  n_tempura_name_mapped <- mapped_name %>% distinct(row_id) %>% nrow()
  n_tempura_any_mapped <- length(unique(c(mapped_taxid$row_id, mapped_name$row_id)))
  n_tempura_unmapped <- nrow(unmapped_final)

  cat(
      "[TEMPURA → GTDB mapping]\n",
      "  TEMPURA rows total:          ", n_tempura_total,        "\n",
      "  With TaxID match (unique):   ", n_tempura_taxid_mapped, "\n",
      "  With name-only match (unique): ", n_tempura_name_mapped,  "\n",
      "  Total matches (unique):     ", n_tempura_any_mapped,   "\n",
      "  With no GTDB match:          ", n_tempura_unmapped,     "\n",
      sep = "")

  list(
    mapped   = mapped_combined,
    unmapped = unmapped_final,
    summary  = list(
      n_tempura_total        = n_tempura_total,
      n_tempura_taxid_mapped = n_tempura_taxid_mapped,
      n_tempura_name_mapped  = n_tempura_name_mapped,
      n_tempura_any_mapped   = n_tempura_any_mapped,
      n_tempura_unmapped     = n_tempura_unmapped
    )
  )
}

#' Wrapper: map full TEMPURA table (Archaea + Bacteria) to GTDB r226
#'   using TaxID first, then names.
#'
#' @param tempura_tbl Full TEMPURA table (both Archaea and Bacteria).
#' @param ar_meta_tbl GTDB Archaea metadata table (ar53_metadata_r226.tsv).
#' @param bac_meta_tbl GTDB Bacteria metadata table (bac120_metadata_r226.tsv).
#'
#' @return List with:
#'   - archaea: list(mapped, unmapped, summary)
#'   - bacteria: list(mapped, unmapped, summary)
#'
map_all_tempura_to_gtdb_r226 <- function(tempura_tbl, ar_meta_tbl, bac_meta_tbl) {

  # Ensure TEMPURA has row_id for tracking
  tempura_with_id <- tempura_tbl %>%
    mutate(row_id = dplyr::row_number())

  tempura_archaea <- tempura_with_id %>%
    filter(superkingdom == "Archaea")

  tempura_bacteria <- tempura_with_id %>%
    filter(superkingdom == "Bacteria")

  cat("\n[Archaea mapping]\n")
  map_arch <- map_tempura_to_gtdb(tempura_archaea, ar_meta_tbl)

  cat("\n[Bacteria mapping]\n")
  map_bac <- map_tempura_to_gtdb(tempura_bacteria, bac_meta_tbl)

  list(
    archaea  = map_arch,
    bacteria = map_bac
  )
}

# Helper: pick the best genome within each TEMPURA row_id
pick_best_gtdb_genome <- function(mapped_tbl) {
  result <- mapped_tbl %>%
    mutate(
      is_gtdb_rep = as.integer(gtdb_representative),  # TRUE=1, FALSE=0
      completeness = coalesce(checkm2_completeness, checkm_completeness),
      contamination = coalesce(checkm2_contamination, checkm_contamination),
      quality_score = coalesce(
        is_gtdb_rep * 1000 + completeness - 5 * contamination,
        is_gtdb_rep * 1000,
        0
      )
    ) %>%
    group_by(row_id) %>%
    slice_max(order_by = quality_score, n = 1, with_ties = FALSE) %>%
    ungroup()
  result
}


