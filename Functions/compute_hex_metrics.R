#==============================================================
# FUNCION: compute_hex_metrics()
# Genera métricas espaciales descriptivas para una malla hexagonal.
# Retorna un dataframe con valores por métrica.
#==============================================================

compute_hex_metrics <- function(
  hex_sf,
  records_col   = "Total_records",
  richness_col  = "Total_sp",
  coverage_col  = "Sampling_coverage",
  majority_prop = 0.8
) {

  if (!inherits(hex_sf, "sf"))
    stop("hex_sf debe ser un objeto sf.")

  #----------------------------------------------------------
  # 1. Área total del estudio
  #----------------------------------------------------------
  area_m2  <- sf::st_area(hex_sf)
  area_ha  <- as.numeric(area_m2) / 10000
  area_tot <- sum(area_ha, na.rm = TRUE)

  recs <- hex_sf[[records_col]]
  rich <- hex_sf[[richness_col]]
  covg <- hex_sf[[coverage_col]]

  recs[is.na(recs)] <- 0
  rich[is.na(rich)] <- NA
  covg[is.na(covg)] <- NA

  #----------------------------------------------------------
  # 2. Área con y sin registros
  #----------------------------------------------------------
  idx_rec <- recs > 0
  area_rec_ha  <- sum(area_ha[idx_rec], na.rm = TRUE)
  pct_area_rec <- 100 * area_rec_ha / area_tot

  #----------------------------------------------------------
  # 3. Celdas que concentran mayoría de registros
  #----------------------------------------------------------
  total_recs <- sum(recs, na.rm = TRUE)

  n_hex_majority   <- NA_integer_
  area_majority_ha <- NA_real_
  pct_area_majority <- NA_real_

  if (total_recs > 0) {
    ord  <- order(recs, decreasing = TRUE)
    recs_sorted <- recs[ord]
    area_sorted <- area_ha[ord]
    cumsum_recs <- cumsum(recs_sorted)

    thr <- majority_prop * total_recs
    idx <- which(cumsum_recs >= thr)[1]

    n_hex_majority   <- idx
    area_majority_ha <- sum(area_sorted[seq_len(idx)])
    pct_area_majority <- 100 * area_majority_ha / area_tot
  }

  #----------------------------------------------------------
  # 4. Estadísticos de riqueza
  #----------------------------------------------------------
  mean_rich   <- mean(rich, na.rm = TRUE)
  median_rich <- stats::median(rich, na.rm = TRUE)

  # outliers por IQR
  outlier_n_cells  <- NA_integer_
  outlier_area_pct <- NA_real_

  rich_clean <- rich[!is.na(rich)]
  if (length(rich_clean) > 0) {
    s <- boxplot.stats(rich_clean)
    out_vals <- s$out
    idx_out  <- rich %in% out_vals & !is.na(rich)
    outlier_n_cells  <- sum(idx_out)
    outlier_area_pct <- 100 * sum(area_ha[idx_out]) / area_tot
  }

  #----------------------------------------------------------
  # 5. Correlación registros – especies
  #----------------------------------------------------------
  idx_cor <- !is.na(rich) & !is.na(recs)
  cor_rec_sp <- NA_real_
  if (sum(idx_cor) > 2 &&
      stats::sd(recs[idx_cor]) > 0 &&
      stats::sd(rich[idx_cor]) > 0) {
    cor_rec_sp <- stats::cor(recs[idx_cor], rich[idx_cor], method = "pearson")
  }

  #----------------------------------------------------------
  # 6. Áreas por cobertura de muestreo
  #----------------------------------------------------------
  idx_cov50 <- covg >= 0.5 & !is.na(covg)
  idx_cov90 <- covg >= 0.9 & !is.na(covg)

  pct_area_cov50 <- 100 * sum(area_ha[idx_cov50]) / area_tot
  pct_area_cov90 <- 100 * sum(area_ha[idx_cov90]) / area_tot

  #----------------------------------------------------------
  # 7. Gini de registros
  #----------------------------------------------------------
  gini_records <- NA_real_
  recs_pos <- recs[recs > 0]

  if (length(recs_pos) > 1) {
    recs_sorted  <- sort(recs_pos)
    n            <- length(recs_sorted)
    gini_records <- (2 * sum(recs_sorted * seq_len(n)) / (n * sum(recs_sorted))) -
      (n + 1) / n
  }

  #----------------------------------------------------------
  # 8. Salida final
  #----------------------------------------------------------
  metrics <- data.frame(
    metric_name = c(
      "total_area_ha",
      "area_with_records_pct",
      "area_with_records_ha",
      "area_without_records_pct",
      "hexes_majority_records_count",
      "hexes_majority_records_area_pct",
      "mean_species_per_hex",
      "median_species_per_hex",
      "outlier_hexes_species_count",
      "outlier_hexes_species_area_pct",
      "records_species_correlation_pearson",
      "area_coverage_ge_0.5_pct",
      "area_coverage_ge_0.9_pct",
      "records_gini_coefficient"
    ),
    metric_description_en = c(
      "Total area of the region in hectares.",
      "Percentage of area with at least one record.",
      "Area (ha) with at least one record.",
      "Percentage of area without records.",
      "Number of hexagons holding most records.",
      "Percentage of total area occupied by those hexagons.",
      "Mean number of species per hexagon.",
      "Median number of species per hexagon.",
      "Number of hexagon richness outliers.",
      "Percentage of total area that these outliers represent.",
      "Pearson correlation between records and species.",
      "Percentage of area with sampling coverage >= 0.5.",
      "Percentage of area with sampling coverage >= 0.9.",
      "Gini coefficient for record spatial inequality."
    ),
    metric_value = c(
      area_tot,
      pct_area_rec,
      area_rec_ha,
      100 - pct_area_rec,
      n_hex_majority,
      pct_area_majority,
      mean_rich,
      median_rich,
      outlier_n_cells,
      outlier_area_pct,
      cor_rec_sp,
      pct_area_cov50,
      pct_area_cov90,
      gini_records
    ),
    stringsAsFactors = FALSE
  )

  return(metrics)
}
