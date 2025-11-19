#==============================================================
# FUNCION: run_SC_abundance_taxon()
# Calcula la cobertura de muestreo (SC) para cada celda
# usando abundancia (frecuencias) por taxón.
# Retorna: un objeto sf con la columna Sampling_coverage.
#==============================================================

run_SC_abundance_taxon <- function(
  taxon_name,
  min_records = 10,
  nboot = 100
) {

  #----------------------------------------------------------
  # 1. Selección de celdas válidas
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
  # 2. Filtrar registros del taxón
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
  # 3. Tabla de salida
  #----------------------------------------------------------
  results_sc <- data.frame(
    id_10km           = integer(),
    Sampling_coverage = numeric()
  )

  pb <- utils::txtProgressBar(min = 0, max = length(id_units), style = 3)
  count <- 0

  #----------------------------------------------------------
  # 4. Iteración por celda
  #----------------------------------------------------------
  for (i in id_units) {

    df_i <- df_taxon %>% dplyr::filter(id_10km == i)

    if (nrow(df_i) < min_records) {
      results_sc <- rbind(
        results_sc,
        data.frame(id_10km = i, Sampling_coverage = NA_real_)
      )
      count <- count + 1
      utils::setTxtProgressBar(pb, count)
      next
    }

    #----------------------------------------------------------
    # 5. Abundancia por taxón (genus)
    #----------------------------------------------------------
    abund_df <- df_i %>%
      dplyr::group_by(genus) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")

    if (nrow(abund_df) == 0) {
      results_sc <- rbind(
        results_sc,
        data.frame(id_10km = i, Sampling_coverage = NA_real_)
      )
      count <- count + 1
      utils::setTxtProgressBar(pb, count)
      next
    }

    abund_vec <- abund_df$n
    names(abund_vec) <- abund_df$genus

    #----------------------------------------------------------
    # 6. iNEXT: cobertura observada
    #----------------------------------------------------------
    sc_obs <- tryCatch({
      r_inext <- iNEXT::iNEXT(
        list(assemblage = abund_vec),
        datatype = "abundance",
        se       = FALSE,
        nboot    = nboot
      )
      as.numeric(r_inext$DataInfo$SC)
    }, error = function(e) NA_real_)

    results_sc <- rbind(
      results_sc,
      data.frame(id_10km = i, Sampling_coverage = sc_obs)
    )

    count <- count + 1
    utils::setTxtProgressBar(pb, count)
  }

  close(pb)

  #----------------------------------------------------------
  # 7. Unir a la capa espacial
  #----------------------------------------------------------
  hex_out <- grid_10km_cut %>%
    dplyr::left_join(results_sc, by = "id_10km")

  return(hex_out)
}
