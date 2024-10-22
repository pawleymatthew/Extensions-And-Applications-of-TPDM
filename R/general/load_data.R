
# EVA Data Challenge ----------------------------------------------------------------

load_eva_c3 <- function(alpha, remove_covariates = TRUE) {
  
  X <- read.csv("data/eva-data-challenge/Coputopia.csv") %>% 
    dplyr::mutate(dplyr::across(Y1:Y3, ~ exp(.x / alpha))) %>%
    dplyr::rename_with(.fn = ~ gsub("Y", "X", .x), .cols = starts_with("Y"))
  
  if (remove_covariates) {
    X <- X %>% dplyr::select(X1, X2, X3) %>% as.matrix()
  }
  
  return(list("X" = X, 
              "n" = nrow(X), 
              "d" = 3,
              "description" = glue::glue("EVA 2023 Data Challenge C3 data on Fréchet({alpha}) margins.")))
}

load_eva_c4 <- function(alpha, order_by_cluster = TRUE, select_clusters = 1:5) {
  
  UtopulaU1 <- read.csv("data/eva-data-challenge/UtopulaU1.csv") %>% as.matrix()
  UtopulaU2 <- read.csv("data/eva-data-challenge/UtopulaU2.csv") %>% as.matrix()
  
  colnames(UtopulaU1) <- paste(colnames(UtopulaU1), "U1", sep = ".")
  colnames(UtopulaU2) <- paste(colnames(UtopulaU2), "U2", sep = ".")
  
  X <- cbind(UtopulaU1, UtopulaU2)
  X <- exp(X / alpha) # standard Gumbel -> Frechet(alpha) margins
  
  clusters <- readRDS("data/eva-data-challenge/Utopula_clusters.RDS")
  clusters <- clusters[clusters %in% select_clusters]
  
  if (order_by_cluster) {
    cluster_order <- clusters %>%
      sort() %>%
      names()
    X <- X[, cluster_order]
  }
  
  return(list("X" = X, 
              "n" = nrow(X), 
              "d" = ncol(X), 
              "clusters" = clusters,
              "description" = glue::glue("EVA 2023 Data Challenge C4 data on Fréchet({alpha}) margins, including components from clusters {glue::glue_collapse(select_clusters,  sep = ', ')}.")))
}

# Kenneth French Finance Data -------------------------------------------------------

# https:// mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library
# date range taken from Meyer and Wintenberger (2023)
# values represent "value-average daily returns of 49 industry portfolios"
# no data transformation applied - see M&W (2023) and Cooley and Thibaud (2019) for related studies

load_kf_49_portfolios <- function(start_date = "1970-01-02", end_date = "2019-12-31") {
  
  X <- "data/kf-finance/49_Industry_Portfolios_Daily.csv" %>%
    read.csv(skip = 9, header = TRUE) %>%
    dplyr::slice_head(n = 25629) %>%
    dplyr::rename(Date = X) %>%
    dplyr::mutate(Date = as.Date(Date, format = "%Y%m%d")) %>%
    dplyr::filter(Date >= start_date & Date <= end_date) %>%
    tibble::column_to_rownames(var = "Date") %>%
    dplyr::mutate(across(.cols = everything(), .fns = function(col) {
      col <- as.numeric(col)
      col <- na_if(col, -99.99)
      return(col)
    })) %>%
    as.matrix()
  
  return(list("X" = X, "n" = nrow(X), "d" = ncol(X)))
}

# Brown-Resnick process ---------------------------------------------

# simulate from a Brown-Resnick process with given variogram across given spatial locations
# default parameters are the simulation study in Supp. Material of Cooley and Thibaud (2019)

sim_cooley_brownresnick <- function(n, alpha = 2, coord = data.frame("x" = rep(1:15, times = 2), "y" = rep(c(0, 15), each = 15)), range = 2.4, smooth = 1.8) {
  
  invisible(capture.output(X <- mvPot::simulBrownResnick(n = n, loc = coord, vario = function(h) (norm(h, "2") / range)^smooth))) # unit Frechet margins
  X <- unlist(X) %>% matrix(nrow = n, ncol = nrow(coord), byrow = TRUE)
  colnames(X) <- paste("X", seq_len(ncol(X)))
  X <- X^(1 / alpha) # Frechet(1) -> Frechet(alpha)
  
  return(list("X" = X, "n" = nrow(X), "d" = ncol(X)))
}

# France rainfall -------------------------------------------------------------------

# https://www.lsce.ipsl.fr/Phocea/Pisp/visu.php?id=109&uid=naveau

load_france_rainfall <- function(alpha = NULL, order_by = bernard_5) {
  
  load("data/france-rainfall/MaxPrecipFallFrance.RData")
  bernard5 <- readRDS(file.path("data", "france-rainfall", "bernard2013-5.rds"))
  bernard7 <- readRDS(file.path("data", "france-rainfall", "bernard2013-7.rds"))
  
  # location info
  coord <- data.frame(
    "id" = colnames(MaxPrecipFallFrance$precip),
    "longitude" = MaxPrecipFallFrance$longitudes,
    "latitude" = MaxPrecipFallFrance$latitudes,
    "bernard_5" = bernard5$clustering,
    "bernard_7" = bernard7$clustering)
  
  # rainfall data
  X <- MaxPrecipFallFrance$precip
  
  # if alpha provided, then standardised to Frechet(alpha) margins
  if (!is.null(alpha)) {
    X <- margins_to_frechet(X, alpha = alpha, type = "empirical")
  }
  
  # order columns of X according to specified variable
  coord <- dplyr::arrange(coord, {{order_by}})
  X <- X[, coord$id]
  
  return(list("X" = X, "coord" = coord, "n" = nrow(X), "d" = ncol(X)))
}

# UK River Flow ---------------------------------------------------------------------

# data emailed to me by Christian Rohrbeck
# river flow data is already transformed to Frechet(alpha=2) margins - see Rohrbeck and Cooley (2023)

load_uk_river_flow <- function(alpha) {
  
  # river flow data
  X <- read.csv(file = "data/uk-river-flow/UK-RiverFlow-Xtilde.csv", header = TRUE) %>% 
    tidyr::drop_na()
  X <- X^(2 / alpha) # Frechet(2) -> Frechet(alpha)
  
  # location info
  coord <- read.csv(file = "data/uk-river-flow/Locations.csv", header = TRUE)
  coord <- data.frame(
    "nrfa_id" = coord$NRFA_ID,
    "name" = coord$Name,
    "tidy_name" = colnames(X),
    "longitude" = coord$Lon, 
    "latitude" = coord$Lat)
  
  return(list("X" = X, "coord" = coord, "n" = nrow(X), "d" = ncol(X)))
}

# Red Sea surface temperature -------------------------------------------------------

# data emailed to me by Christian Rohrbeck
# temperature values have been preprocessed so that margins are stationary (ask Christian what the distribution is?)

load_red_sea_temp <- function(alpha) {

  coord <- readRDS("data/red-sea-temperature/Locations_Subset.Rda") %>%
    rename(longitude = lon, latitude = lat) %>%
    mutate(region = case_when(
      longitude > 34.5 & longitude < 38.5 & latitude > 19.5  ~ "north",
      longitude > 39 & longitude < 42 & latitude < 20  ~ "south",
      .default = "none"))
  
  # temperature data - weekly maxima
  sst.gauss <- readRDS("data/red-sea-temperature/StatMarginsRedSeaAnomaly.RDS")
  X <- sst.gauss %>% 
    mutate(week = (row_number() - 1) %/% 7) %>%
    group_by(week) %>%
    summarise_all(max) %>%
    select(-week) %>%
    as.matrix() %>%
    margins_to_frechet(alpha = alpha, type = "empirical")

  return(list("X" = X, "coord" = coord, "n" = nrow(X), "d" = ncol(X)))
}


# Ecoli -----------------------------------------------------------------------------

load_ecoli <- function() {
  read.table("data/ecoli/ecoli.data", dec = ",", header = FALSE) %>%
    as_tibble() %>%
    set_colnames(c("sequence", "mcg", "gvh", "lip", "chg", "aac", "alm1", "alm2", "class")) %>%
    mutate(across(mcg:alm2, .fns = as.numeric)) %>%
    mutate(class = case_when(class == "im" ~ "1", .default = "2")) %>% # im = 1, not im = 2
    select(-lip, -chg) %>%
    group_by(class) %>%
    group_split() %>%
    lapply(function(df) {
      df[, c("mcg", "gvh", "aac", "alm1", "alm2")] <- margins_to_frechet(X = as.matrix(df[, c("mcg", "gvh", "aac", "alm1", "alm2")]), alpha = 1, type = "empirical")
      return(df)
    }) %>%
    bind_rows() %>%
    mutate(class = as.factor(class))
}


# NFL Combine data ------------------------------------------------------------------

load_nfl_combine <- function(standardise = "byclass", frechet = "byclass", alpha = 1) {
  
  tmp <- list.files("data/nfl-combine/", pattern="\\.csv$", full.names = TRUE) %>%
    lapply(read_csv, show_col_types = FALSE) %>%
    bind_rows(.id = "year") %>%
    mutate(year = 1999 + as.integer(year)) %>%
    select(-School) %>%
    rowwise() %>%
    mutate(Ht = sum(c(0.3048, 0.0254) * as.numeric(unlist(str_split(Ht, pattern = "-"))))) %>% # height in metres
    ungroup() %>%
    drop_na() %>%
    rename(FortyYard = `40yd`, BroadJump = `Broad Jump`, ThreeCone = `3Cone`) %>%
    mutate(FortyYard = 40 / FortyYard, # total time [sec] -> average speed [yards/sec]
           ThreeCone = 30 / ThreeCone, 
           Shuttle = 20 / Shuttle) %>%
    mutate(class = case_when(
      Pos %in% c("C", "OG", "OT", "DE", "DT") ~ "1", #  1 = on-the-line, 2 = off-the-line
      .default = "2")) %>%
    mutate(class = as.factor(class))

  if (!isFALSE(standardise)) {
    if (standardise == "byclass") {
      tmp <- tmp %>%
        group_by(class) %>%
        group_split() %>%
        lapply(function(df) {
          df %>% mutate(across(.cols = Ht:Shuttle,
                               .fns = ~ scale(.x, center = TRUE, scale = TRUE)))
        }) %>%
        bind_rows()
    } else {
      tmp <- tmp %>% 
        mutate(across(.cols = Ht:Shuttle,
                      .fns = ~ scale(.x, center = TRUE, scale = TRUE)))
    }
  }
  
  if (!isFALSE(frechet)) {
    if (frechet == "byclass") {
      tmp <- tmp %>%
        group_by(class) %>%
        group_split() %>%
        lapply(function(df) {
          df %>% mutate(across(.cols = Ht:Shuttle,
                               .fns = ~ margins_to_frechet(X = .x, alpha = alpha, type = "empirical")))
        }) %>%
        bind_rows()
    } else {
      tmp <- tmp %>% 
        mutate(across(.cols = Ht:Shuttle,
                      .fns = ~ margins_to_frechet(X = .x, alpha = alpha, type = "empirical")))
    }
  }
  
  tmp <- mutate(tmp, across(.cols = Ht:Shuttle, .fns = ~ as.numeric(.x)))

  # if (standardise) {
  #   tmp <- tmp %>%
  #     group_by(class) %>%
  #     group_split() %>%
  #     lapply(function(df) {
  #       df <- df %>% mutate(across(.cols = Ht:Shuttle,
  #                                  .fns = ~ scale(.x, center = TRUE, scale = TRUE)))
  #       if (frechet) {
  #         df <- df %>% mutate(across(.cols = Ht:Shuttle,
  #                                    .fns = ~ margins_to_frechet(X = .x, alpha = 1, type = "empirical")))
  #       }
  #       return(df)
  #     }) %>%
  #     bind_rows() %>%
  #     mutate(across(.cols = Ht:Shuttle, .fns = ~ as.numeric(.x)))
  # }
  return(tmp)
}

