# analysis/functions/tempura_helpers.R
# Helper functions for working with the TEMPURA growth-temperature database.
# TEMPURA reference: Sato et al. 2020, Database of Growth TEMPeratures of Usual and RAre Prokaryotes. 
# doi: 10.1264/jsme2.ME20074

#' Load TEMPURA and filter to one superkingdom
#'
#' @param tempura_path Character path to TEMPURA CSV file.
#' @param superkingdom_name Character, e.g. "Archaea" or "Bacteria".
#' @return A list with:
#'   - data_focal: TEMPURA rows for the chosen superkingdom
#'   - with_accession: subset with non-NA assembly_or_accession
#'   - without_accession: subset lacking assembly_or_accession
#'   - taxid_missing: rows in without_accession with missing taxonomy_id
#'   - summary: named list of counts
#'
load_tempura_superkingdom <- function(tempura_path, superkingdom_name) {
  # Define EXACT columns you want
  tempura_cols <- c("genus_and_species", "taxonomy_id", "strain", "superkingdom", 
                    "phylum", "class", "order", "family", "genus", 
                    "assembly_or_accession", "Tmin", "Topt_ave", "Topt_low", 
                    "Topt_high", "Tmax", "Tmax_Tmin")
  
  # Read + subset to ONLY needed columns
  tempura <- readr::read_csv(tempura_path) %>%
    dplyr::select(all_of(tempura_cols))
  
  # Filter + add row_id
  data_focal <- tempura %>%
    dplyr::filter(superkingdom == superkingdom_name) %>%
    dplyr::mutate(row_id = dplyr::row_number())
  
  # Split by accession presence (QC)
  with_accession <- data_focal %>% dplyr::filter(!is.na(assembly_or_accession))
  without_accession <- data_focal %>% dplyr::filter(is.na(assembly_or_accession))
  taxid_missing <- without_accession %>% dplyr::filter(is.na(taxonomy_id))
  
  summary_counts <- list(
    superkingdom      = superkingdom_name,
    n_total           = nrow(data_focal),
    n_with_accession  = nrow(with_accession),
    n_without_access  = nrow(without_accession),
    n_taxid_missing   = nrow(taxid_missing)
  )
  
  cat("[TEMPURA CLEAN]", superkingdom_name, nrow(data_focal), "rows\n")
  
  list(
    data_focal        = data_focal,        # MAIN OUTPUT (use this!)
    with_accession    = with_accession,    # QC splits
    without_accession = without_accession,
    taxid_missing     = taxid_missing,
    summary           = summary_counts
  )
}

#' Load COMPLETE TEMPURA (Archaea + Bacteria combined)
load_tempura_complete <- function(tempura_path) {
  ar_tempura <- load_tempura_superkingdom(tempura_path, "Archaea")$data_focal
  bac_tempura <- load_tempura_superkingdom(tempura_path, "Bacteria")$data_focal
  bind_rows(ar_tempura, bac_tempura)
}