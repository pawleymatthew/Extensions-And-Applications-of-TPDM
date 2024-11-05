# Load required package
if (!require("ggtern")) install.packages("ggtern")
library(ggtern)

# Calculate Aitchison distance for each grid point
calculate_aitchison_distance <- function(point, ref) {
  clr_point <- log(point) - mean(log(point))
  clr_ref <- log(ref) - mean(log(ref))
  dist <- sqrt(sum((clr_point - clr_ref)^2))
  return(dist)
}

calculate_euclidean_distance <- function(point, ref) {
  dist <- sqrt(sum((point - ref)^2))
  return(dist)
}

ref_list <- list("centre" = c(1/3, 1/3, 1/3), "off_centre" = c(0.105, 0.4475, 0.4475))

res <- lapply(ref_list, function(ref) {
  
  # Generate a grid of points within the simplex
  grid_points <- expand_grid(A = seq(0.005, 0.995, length = 99),
                             B = seq(0.005, 0.995, length = 99)) %>%
    filter(A + B <= 1) %>%
    mutate(C = 1 - A - B) %>%
    filter(C >= 0.005) %>%
    rowwise() %>%
    mutate(point = list(c(A, B, C))) %>%
    rowwise() %>%
    mutate(aitchison_distance = calculate_aitchison_distance(point, ref)) %>%
    mutate(euclidean_distance = calculate_euclidean_distance(point, ref)) %>%
    ungroup()
  
  return(grid_points)
  
}) %>% 
  bind_rows(.id = "ref_point") %>%
  mutate(ref_point_coord = list(ref_list[[ref_point]])) %>%
  pivot_longer(cols = ends_with("distance"), names_to = "distance_type", values_to = "distance")

saveRDS(res, file.path("scripts", "compositional", "results", "aitchison_distance.RDS"))

