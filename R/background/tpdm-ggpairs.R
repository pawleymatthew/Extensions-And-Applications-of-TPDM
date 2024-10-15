library(GGally)
library(scales)

tpdm_max_linear_ggpairs <- function(data, Sigma, V, k) {
  
  custom_scatter <- function(data, mapping, ...) {
    
    vecu_sigma <- set_names(vecu(Sigma), colnames(V))
    bivnorm_mean <- vecu_sigma[c(quo_name(mapping$x), quo_name(mapping$y))]
    bivnorm_cov <- V[c(quo_name(mapping$x), quo_name(mapping$y)),
                     c(quo_name(mapping$x), quo_name(mapping$y))]
    bivnorm_cov_eigen <- eigen(bivnorm_cov / k)
    ellipse_a <- sqrt(qchisq(0.95, 2) * bivnorm_cov_eigen$values[1])
    ellipse_b <- sqrt(qchisq(0.95, 2) * bivnorm_cov_eigen$values[2])
    ellipse_angle <- atan2(bivnorm_cov_eigen$vectors[2, 1], bivnorm_cov_eigen$vectors[1, 1])
    
    p <- ggplot(data = data, mapping = mapping) + 
      geom_point(colour = "grey", size = 0.3) +
      annotate(geom = "point", x = bivnorm_mean[1], y = bivnorm_mean[2], shape = 4, colour = "blue") +
      annotate(geom = "point", x = mean(eval_data_col(data, mapping$x)), y = mean(eval_data_col(data, mapping$y)), shape = 4, colour = "red") +
      stat_ellipse(colour = "red") +
      geom_ellipse(aes(x0 = bivnorm_mean[1], y0 = bivnorm_mean[2], 
                       a = ellipse_a, b = ellipse_b, angle = ellipse_angle),
                   colour = "blue", size = 0.1) +
      scale_x_continuous(breaks = breaks_extended(n = 3)) +
      scale_y_continuous(breaks = breaks_extended(n = 3)) +
      theme_light()
    p
  }
  custom_hist <- function(data, mapping, ...){
    
    vecu_sigma <- set_names(vecu(Sigma), colnames(V))
    dnorm_mean <- vecu_sigma[quo_name(mapping$x)]
    dnorm_sd <- sqrt(V[quo_name(mapping$x), quo_name(mapping$x)] / k)
    
    p <- ggplot(data = data, mapping = mapping) + 
      geom_histogram(aes(y = after_stat(density)), fill = "red", alpha = 0.4) +
      geom_vline(xintercept = dnorm_mean, colour = "blue", linetype = "dashed") +
      geom_function(fun = dnorm, args = list(mean = dnorm_mean, sd = dnorm_sd), colour = "blue") +
      scale_x_continuous(breaks = breaks_extended(n = 3)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
      theme_light()
    p
  }
  custom_corr <- function(data, mapping, ...) {
    
    cov_emp <- (k * cov(eval_data_col(data, mapping$x), eval_data_col(data, mapping$y))) %>%
      scientific(digits = 2, trim = FALSE) %>%
      gsub("e", " %*% 10^{", .) %>%
      paste0("}")
    
    cov_true <- V[quo_name(mapping$x), quo_name(mapping$y)] %>%
      scientific(digits = 2, trim = FALSE) %>%
      gsub("e", " %*% 10^{", .) %>%
      paste0("}")
    
    p <- ggplot(data = data, mapping = mapping) +
      theme_void() +
      annotate("text", x = mean(eval_data_col(data, mapping$x)), y = mean(eval_data_col(data, mapping$y)), 
               hjust = 0.5, vjust = 1, 
               label = cov_emp, parse = TRUE, 
               color = "red") +
      annotate("text", x = mean(eval_data_col(data, mapping$x)), y = mean(eval_data_col(data, mapping$y)), 
               hjust = 0.5, vjust = 0, 
               label = cov_true, parse = TRUE, 
               color = "blue")
    return(p)
  }
  
  p <- ggpairs(data, 
               diag = list(continuous = wrap(custom_hist, Sigma = Sigma, V = V)),
               lower = list(continuous = wrap(custom_scatter, Sigma = Sigma, V = V)),
               upper = list(continuous = wrap(custom_corr, k = 100)),
               labeller = label_parsed,
               progress = FALSE)
  
  return(p)
}