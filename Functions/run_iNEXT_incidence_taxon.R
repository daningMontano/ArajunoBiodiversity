#==============================================================
# FUNCION: run_iNEXT_incidence_taxon()
# Calcula métricas de iNEXT (incidencia) para un taxón específico.
# Retorna:
#   - result: objeto iNEXT completo
#   - metrics: dataframe con cobertura, riqueza observada,
#              riqueza extrapolada y unidades observadas/extrapoladas
#==============================================================

run_iNEXT_incidence_taxon <- function(
  taxon_name,
  min_records = 10,
  nboot = 100,
  save_path = NULL
) {

  #----------------------------------------------------------
  # 1. Filtrar celdas válidas según el taxón
  #----------------------------------------------------------
  id_units <- results_taxa_total %>%
    dplyr::filter(
      Group == taxon_name,
      Total_records >= min_records,
      !is.na(id_10km)
    ) %>%
    dplyr::pull(id_10km) %>%
    unique()

  if (length(id_units) == 0)
    stop("No hay celdas que cumplan el filtro")

  cat("Celdas filtradas:", length(id_units), "\n")

  #----------------------------------------------------------
  # 2. Filtrar registros en esas celdas
  #----------------------------------------------------------
  df_taxon <- recors_grid_10km %>%
    dplyr::filter(
      group == taxon_name,
      id_10km %in% id_units
    )

  if (nrow(df_taxon) == 0)
    stop("No hay registros del taxón en estas celdas")

  cat("Registros encontrados:", nrow(df_taxon), "\n")

  #----------------------------------------------------------
  # 3. Matriz de incidencia (especies x celdas)
  #----------------------------------------------------------
  incid_wide <- df_taxon %>%
    dplyr::distinct(id_10km, scientificName) %>%
    dplyr::mutate(val = 1L) %>%
    tidyr::pivot_wider(
      id_cols    = id_10km,
      names_from = scientificName,
      values_from = val,
      values_fill = 0
    )

  ids <- incid_wide$id_10km

  incid_mat <- incid_wide %>%
    dplyr::select(-id_10km) %>%
    as.data.frame()

  # transponer: especies como filas
  mat_final <- as.data.frame(t(incid_mat))

  cat("Especies:", nrow(mat_final),
      "| Celdas:", ncol(mat_final), "\n")

  #----------------------------------------------------------
  # 4. Ejecutar iNEXT
  #----------------------------------------------------------
  result_inext <- iNEXT::iNEXT(
    list(assemblage = mat_final),
    q = 0,
    datatype = "incidence_raw",
    se = TRUE,
    conf = 0.95,
    nboot = nboot
  )

  if (!is.null(save_path)) {
    saveRDS(result_inext, save_path)
  }

  #----------------------------------------------------------
  # 5. Extraer métricas clave
  #----------------------------------------------------------

  sc <- result_inext$DataInfo$SC

  sp_obs <- result_inext$iNextEst$size_based %>%
    dplyr::filter(Method == "Observed") %>%
    dplyr::pull(qD)

  sp_ext <- result_inext$iNextEst$size_based %>%
    dplyr::filter(Method == "Extrapolation") %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::pull(qD)

  units_obs <- result_inext$iNextEst$size_based %>%
    dplyr::filter(Method == "Observed") %>%
    dplyr::pull(t)

  units_ext <- result_inext$iNextEst$size_based %>%
    dplyr::filter(Method == "Extrapolation") %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::pull(t)

  metrics <- data.frame(
    sample_coverage_observed = sc,
    species_observed         = sp_obs,
    species_extrapolated     = sp_ext,
    units_observed           = units_obs,
    units_extrapolated       = units_ext
  )

  cat("SC:", sc,
      "| Sp_obs:", sp_obs,
      "| Units_obs:", units_obs,
      "| Units_ext:", units_ext,
      "| Sp_ext:", sp_ext, "\n")

  #----------------------------------------------------------
  # 6. Retornar
  #----------------------------------------------------------
  return(list(
    result  = result_inext,
    metrics = metrics
  ))
}
