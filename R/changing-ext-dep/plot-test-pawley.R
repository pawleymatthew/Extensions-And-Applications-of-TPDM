plot_test_pawley <- function(data, variable_scheme = "colour", format = c(2, 2), legend_position = "none") {
  
  n_blocks <- max(data$block)
  D <- data$variables %>% unique() %>% length()
  cm_cv <- quantile(bb_L2, probs = 0.95^(2/D)) %>% as.numeric()
  ks_cv <- CPAT:::qkolmogorov(0.95^(2/D))
  
  if (D > 6 & variable_scheme == "facet") {
    print("variable_scheme = facet is only available for d<=4.")
    variable_scheme <- "colour"
  }
  
  # add a block zero row to facilitate plotting
  data <- data %>%
    rbind(expand_grid(block = 0, sigma = 0, psi = 0, z = 0, variables = unique(pull(., variables)), ks = 0, cm = 0))
  
  if (variable_scheme == "facet") { # use facet_grid
    
    p_sigma <- filter(data, block > 0) %>%
      rbind(tail(., n = D) %>% mutate(block = block + 1)) %>%
      ggplot(aes(x = (block - 1) / n_blocks, y = sigma)) +
      facet_grid(~ variables)
    
    p_psi <- data %>%
      ggplot(aes(x = block / n_blocks, y = psi)) +
      facet_grid(~ variables) 
    
    p_z <- data %>%
      ggplot(aes(x = block / n_blocks, y = z)) +
      facet_grid(~ variables) 
    
    p_cv <- data %>%
      pivot_longer(cols = c(ks, cm), names_to = "type", values_to = "statistic") %>%
      mutate(cv = case_when(type == "cm" ~ cm_cv, .default = ks_cv)) %>%
      ggplot(aes(x = block / n_blocks, y = statistic)) +
      facet_grid(type ~ variables, scales = "free", labeller = labeller(.default = toupper))
    
  } else if (variable_scheme == "colour") { # use colours
    
    p_sigma <- filter(data, block > 0) %>%
      rbind(tail(., n = D) %>% mutate(block = block + 1)) %>%
      ggplot(aes(x = (block - 1) / n_blocks, y = sigma, colour = variables)) +
      scale_colour_discrete_sequential("ag_GrnYl")
    
    p_psi <- data %>%
      ggplot(aes(x = block / n_blocks, y = psi, colour = variables)) +
      scale_colour_discrete_sequential("ag_GrnYl")
    
    p_z <- data %>%
      ggplot(aes(x = block / n_blocks, y = z, colour = variables)) +
      scale_colour_discrete_sequential("ag_GrnYl")
    
    p_cv <- data %>%
      pivot_longer(cols = c(ks, cm), names_to = "type", values_to = "statistic") %>%
      mutate(cv = case_when(type == "cm" ~ cm_cv, .default = ks_cv)) %>%
      ggplot(aes(x = block / n_blocks, y = statistic, colour = variables)) +
      facet_grid(type ~ ., scales = "free", labeller = labeller(.default = toupper)) +
      scale_colour_discrete_sequential("ag_GrnYl")
    
  } else {
    stop("Specify whether variable pairs should be represented by colour or facet.")
  }
  
  p_sigma <- p_sigma +
    geom_step() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_extended(n = 4)) +
    xlab("Time") +
    ylab(expression(hat(sigma)[ij](t))) +
    labs(colour = "Variables") +
    theme_light()
  
  p_psi <- p_psi +
    geom_line() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_extended(n = 4)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02), add = c(0, 0))) +
    xlab("Time") +
    ylab(expression(hat(psi)[ij](t))) +
    labs(colour = "Variables") +
    theme_light()
  
  p_z <- p_z +
    geom_line() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_extended(n = 4)) +
    xlab("Time") +
    ylab(expression(hat(Z)[ij](t))) +
    labs(colour = "Variables") +
    theme_light()
  
  p_cv <- p_cv +
    geom_line() +
    geom_hline(aes(yintercept = cv), colour = "black", linetype = "dashed") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_extended(n = 4)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1), add = c(0, 0)), breaks = breaks_extended(n = 4)) +
    xlab("Time") +
    ylab("Test statistic process") +
    labs(colour = "Variables") +
    theme_light()
  
  p <- ggarrange(p_sigma, p_psi, p_z, p_cv, nrow = format[1], ncol = format[2], common.legend = TRUE, legend = legend_position)
  return(p)
  
}

