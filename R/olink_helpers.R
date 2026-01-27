#' Calculate Coefficient of Variation (CV) for Olink Data
#'
#' Calculates intra-plate and inter-plate CVs based on Sample Controls.
#' Note: Uses the logarithmic transformation method for CV calculation 
#' characteristic of Olink data.
#'
#' @param dat A data frame containing Olink NPX data.
#' @return A data frame with joined intra and inter-plate CV metrics.
#' @export
calculate_olink_cv <- function(dat) {
  # 1. Filter for controls and calculate intra-plate metrics
  sample_controls <- dat %>% 
    dplyr::filter(SampleType == "SAMPLE_CONTROL") %>% 
    dplyr::select(OlinkID, PlateID, NPX)
  
  intra_cv <- sample_controls |> 
    dplyr::group_by(OlinkID, PlateID) |> 
    dplyr::summarise(
      NPX_median = mean(NPX, na.rm = TRUE),
      MyIntraCV = 100 * sqrt(exp((log(2) * sd(NPX, na.rm = TRUE))^2) - 1),
      .groups = 'drop'
    ) 
  
  # 2. Calculate inter-plate metrics
  inter_cv <- intra_cv |> 
    dplyr::group_by(OlinkID) |> 
    dplyr::summarise(
      MyInterCV = 100 * sqrt(exp((log(2) * sd(NPX_median, na.rm = TRUE))^2) - 1),
      MyInterCV2 = mean(NPX_median, na.rm = TRUE) / sd(NPX_median, na.rm = TRUE),
      .groups = 'drop' 
    )
  
  # 3. Join back to original data
  rv <- dat |> 
    dplyr::inner_join(intra_cv |> dplyr::select(-NPX_median), by = c("OlinkID", "PlateID")) |> 
    dplyr::inner_join(inter_cv, by = "OlinkID")
  
  return(rv)
}

#' Plot Histogram of Intra-plate CVs
#'
#' @param dat Data frame containing IntraCV and LOD columns.
#' @param out_fn Character. Output file path for the PDF.
#' @param width Numeric. Plot width.
#' @param height Numeric. Plot height.
#' @export
plot_intra_cv_histogram <- function(dat, out_fn, width = 10, height = 8) {
  require(ggplot2)
  
  p <- dat |> 
    dplyr::filter(!is.na(IntraCV)) |> 
    ggplot(aes(x = IntraCV)) +
    geom_histogram(
      binwidth = 0.5, 
      fill = "#0072B2", 
      alpha = 0.8
    ) + 
    facet_grid(Block ~ (NPX > LOD)) +
    theme_minimal() +
    labs(title = "Intra-CV Distribution", x = "Intra-CV (%)", y = "Count")
  
  ggsave(out_fn, p, width = width, height = height)
}

#' Check NPX Data Integrity
#'
#' Performs various checks on Olink NPX data frames including 
#' OID conformity, missing values, and QC warnings.
#'
#' @param df A data frame containing Olink NPX data.
#' @return A list of check results (non-conforming OIDs, NAs, duplicates).
#' @export
check_npx_integrity <- function(df) {
  df_colnames <- colnames(df)
  
  # Check data type
  if ("NPX" %in% df_colnames) {
    data_type <- "NPX"
  } else if ("Quantified_value" %in% df_colnames) {
    data_type <- "Quantified_value"
  } else {
    stop("Neither NPX nor Quantified_value column present in data.")
  }
  
  # OlinkID Check
  if (!("OlinkID" %in% df_colnames)) stop("OlinkID column missing.")
  
  non_conforming_oid <- df |> 
    dplyr::distinct(OlinkID) |> 
    dplyr::filter(!grepl("OID[0-9]{5}", OlinkID)) |> 
    dplyr::pull(OlinkID)
  
  # Uniprot ID Checks
  uniprot_replacements <- find_uniprot_mismatches(df)
  
  # Duplicate SampleID/OlinkID pairs
  duplicate_rows <- df |> 
    dplyr::select(SampleID, OlinkID) |> 
    duplicated()
  
  duplicate_samples <- if (any(duplicate_rows)) unique(df$SampleID[duplicate_rows]) else character(0)
  
  # Summarize results
  return(list(
    data_type = data_type,
    non_conforming_oid = non_conforming_oid,
    duplicate_samples = duplicate_samples,
    uniprot_replacements = uniprot_replacements
  ))
}

#' Identify Mismatched UniProt IDs for OlinkIDs
#'
#' @param df Olink data frame.
#' @return Data frame of OlinkIDs with inconsistent UniProt mappings.
find_uniprot_mismatches <- function(df) {
  assay_ids <- df |> 
    dplyr::select(dplyr::any_of(c("OlinkID", "UniProt"))) |> 
    dplyr::distinct()
  
  duplicated_oids <- assay_ids$OlinkID[duplicated(assay_ids$OlinkID)]
  
  if (length(duplicated_oids) > 0) {
    # Logic to identify which UniProt to keep (usually the first)
    # Returns mapping for replacement
    cli::cli_warn("{length(unique(duplicated_oids))} OlinkID(s) have multiple UniProt IDs.")
  }
  return(duplicated_oids)
}

#' Process Olink Data for PCA
#'
#' Filters low variance assays, pivots to wide format, and imputes missing values with the median.
#'
#' @param dat Olink NPX data.
#' @param sam Sample metadata.
#' @return A list containing the wide 'dat' matrix and aligned 'sam' metadata.
#' @export
prepare_pca_data <- function(dat, sam) {
  # Remove assays with zero variance
  df_filtered <- dat |> 
    dplyr::group_by(OlinkID, PlateID) |> 
    dplyr::mutate(assay_var = var(NPX, na.rm = TRUE)) |> 
    dplyr::ungroup() |> 
    dplyr::filter(!(assay_var == 0 | is.na(assay_var))) |> 
    dplyr::select(-assay_var)
  
  # Pivot to wide format
  df_wide <- df_filtered |> 
    dplyr::select(SampleID, OlinkID, NPX) |> 
    dplyr::filter(!is.na(NPX)) |> 
    tidyr::pivot_wider(names_from = OlinkID, values_from = NPX, names_sort = TRUE) |> 
    tibble::column_to_rownames("SampleID") |> 
    dplyr::mutate(dplyr::across(where(is.numeric), ~tidyr::replace_na(., median(., na.rm = TRUE))))
  
  # Align metadata
  sam_aligned <- sam |> 
    dplyr::mutate(ID_Tmp = SampleID) |> 
    tibble::column_to_rownames("ID_Tmp")
  
  return(list(dat = df_wide, sam = sam_aligned[rownames(df_wide), ]))
}

#' Generate Olink PCA Plot
#'
#' Wraps processing, PCA calculation, and plotting into one call.
#'
#' @param dat NPX data.
#' @param sam Metadata.
#' @param color_by Variable to color points by.
#' @param facet_by Optional variable to facet the plot.
#' @export
run_olink_pca_workflow <- function(dat, sam, color_by = "SampleType", facet_by = "PlateID", ncol=2) {
  # 1. Prepare
  prepared <- prepare_pca_data(dat, sam)
  
  # 2. Compute PCA
  pca_res <- prcomp(prepared$dat, scale. = TRUE, center = TRUE)
  pca_vars <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  
  # 3. Join with metadata
  plot_df <- cbind(prepared$sam, pca_res$x[, 1:5])
  
  # 4. Plot
  p <- ggplot2::ggplot(plot_df, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", pca_vars[1]),
      y = sprintf("PC2 (%.1f%%)", pca_vars[2])
    ) +
    geom_vline(xintercept = 0, lty=2, alpha=0.3) +  
    geom_hline(yintercept = 0, lty=2, alpha=0.3) + 
    coord_fixed() +
    ggsci::scale_color_jco() +
    theme(plot.margin=unit(c(1,1,1,1),"line"), panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank()) # the theem to have a clear background and border enclose the panel
  
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_by)), ncol=ncol)
  }
  
  return(p)
}

#' Normalize Olink Data (Intensity or Plate Control)
#'
#' @param dat Olink data with ExtNPX values.
#' @return Data frame with custom normalized NPX columns.
#' @export
normalize_olink_custom <- function(dat) {
  dat |> 
    dplyr::group_by(PlateID, OlinkID) |> 
    dplyr::mutate(
      pc_median = median(ExtNPX[SampleType == "PLATE_CONTROL"], na.rm = TRUE),
      sample_median = median(ExtNPX[SampleType == "SAMPLE"], na.rm = TRUE),
      NPX_PC_Norm = ifelse(SampleQC %in% c("FAIL", "NA"), NA, ExtNPX - pc_median),
      NPX_Intensity_Norm = ifelse(SampleQC %in% c("FAIL", "NA"), NA, ExtNPX - sample_median)
    ) |> 
    dplyr::ungroup()
}

#' Recode NPX Values Based on LOD and Missingness
#'
#' Sets NPX to NA if below LOD or if >50% of the assay/plate is missing.
#'
#' @param dat Olink data.
#' @export
recode_low_quality_npx <- function(dat) {
  dat |> 
    dplyr::mutate(NPX = ifelse(NPX < LOD, NA, NPX)) |> 
    dplyr::group_by(PlateID, OlinkID) |> 
    dplyr::mutate(
      too_many_nas = mean(is.na(NPX)) > 0.5,
      NPX = ifelse(too_many_nas, NA, NPX)
    ) |> 
    dplyr::ungroup() |> 
    dplyr::select(-too_many_nas) |> 
    dplyr::group_by(OlinkID) |> 
    dplyr::mutate(
      global_too_many_nas = mean(is.na(NPX)) > 0.5,
      NPX = ifelse(global_too_many_nas, NA, NPX)
    ) |> 
    dplyr::filter(sum(!is.na(NPX)) > 0) |> 
    dplyr::ungroup() |> 
    dplyr::select(-global_too_many_nas)
}

#' Inspect Olink QC Flags for Assays and Samples
#' 
#' This function identifies assays and samples that do not have a "PASS" QC status,
#' providing a detailed breakdown of SampleQC status across different blocks.
#'
#' @param dat A tibble containing Olink NPX data.
#' @return A list containing dataframes for flagged assays, sample IDs, and a block-wise sample QC matrix.
inspect_olink_qc <- function(dat) {
  
  # 1. Identify Assays that are not "PASS" (excluding controls)
  # Filtering for AssayType=="assay" ensures we focus on actual targets
  flagged_assays <- dat |> 
    filter(AssayType == "assay") |> 
    select(OlinkID, Assay, Block, AssayQC) |> 
    distinct() |> 
    filter(AssayQC != "PASS")
  
  # 2. Identify unique SampleIDs that are not "PASS" or "NA"
  # Based on Olink definitions, "NA" usually refers to excluded assays, not failed samples
  flagged_sample_ids <- dat |> 
    filter(!SampleQC %in% c("NA", "PASS")) |> 
    pull(SampleID) |> 
    unique()
  
  # 3. Create a block-wise matrix of SampleQC status for flagged samples
  # This helps visualize if a sample failed across all blocks or just one
  flagged_samples_matrix <- dat |> 
    filter(SampleID %in% flagged_sample_ids) |>
    group_by(SampleID, Block, SampleQC) |> 
    summarise(N = n(), .groups = 'drop') |>
    mutate(status_label = paste0(SampleQC, " (n=", N, ")")) |>
    group_by(SampleID, Block) |> 
    summarise(QC = stringr::str_flatten(status_label, collapse = ", "), .groups = 'drop') |> 
    tidyr::pivot_wider(names_from = Block, values_from = QC)
  
  # Return results as a named list
  list(
    flagged_assays = flagged_assays, 
    flagged_sample_ids = flagged_sample_ids,
    flagged_samples_matrix = flagged_samples_matrix
  )
}


#' Calculate Intra-CV using Log-Normal Transformation
#' 
#' @param dat Data frame containing Olink NPX results.
#' @param cgr_sam Mapping file for biological samples.
#' @param olink_sam Data frame containing control sample definitions.
#' @param npx_col Bare name of the NPX column.
#' @param lod_col Bare name of the LOD column.
#' @param treat_LOD_as_NA Logical; if TRUE, values < LOD are ignored.
calc_intraCV_nse <- function(dat, 
                             cgr_sam, 
                             olink_sam, 
                             npx_col = NPX, 
                             lod_col = LOD,
                             treat_LOD_as_NA = FALSE) {
  
  # --- 1. Prepare Subject/Sample Mapping ---
  # Combine biological sample mapping with control sample definitions
  sam_subj <- cgr_sam |> 
    dplyr::select(SampleID, subsample_id) |> 
    dplyr::bind_rows(
      olink_sam |> 
        dplyr::filter(grepl("_CONTROL", SampleType)) |> 
        dplyr::select(SampleID, subsample_id = SampleType)
    ) |> 
    dplyr::distinct() # Prevent duplicate joins
  
  # --- 2. Data Preparation and LOD Handling ---
  dat_processed <- dat |>
    dplyr::mutate(
      NPX_temp = {{npx_col}}, 
      LOD_temp = {{lod_col}}
    )
  
  if (treat_LOD_as_NA) {
    dat_processed <- dat_processed |>
      dplyr::mutate(NPX_temp = dplyr::if_else(NPX_temp < LOD_temp, NA_real_, NPX_temp))
  }
  
  # --- 3. Calculate IntraCV ---
  # Calculation uses the log-normal CV formula: 100 * sqrt(exp(sigma^2) - 1)
  # Since NPX is log2, we convert sd to natural log scale: sd * log(2)
  intraCV <- dat_processed |> 
    dplyr::inner_join(sam_subj, by = "SampleID") |>
    dplyr::group_by(OlinkID, Block, AssayType, PlateID, subsample_id) |>
    # Ensure at least 2 replicates exist for a valid SD calculation
    dplyr::filter(sum(!is.na(NPX_temp)) > 1) |>
    dplyr::summarise(
      IntraCV = 100 * sqrt(exp((log(2) * sd(NPX_temp, na.rm = TRUE))^2) - 1),
      n_replicates = n(),
      .groups = 'drop'
    )
  
  return(intraCV)
}