plot_test_drees <- function(data, format = c(1, 3)) {
  
  n_blocks <- max(data$block)
  n_sets <- data$set %>% unique() %>% length()
  
  # add a block zero row to facilitate plotting
  data <- data %>%
    rbind(expand_grid(block = 0, S = 0, IS = 0, z = 0, set = unique(pull(., set)), ks = 0, cm = 0))

  p_S <- filter(data, block > 0) %>%
    rbind(tail(., n = n_sets) %>% mutate(block = block + 1)) %>%
    ggplot(aes(x = (block - 1) / n_blocks, y = S, colour = set)) +
    scale_colour_discrete_sequential("ag_GrnYl") +
    geom_step() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_pretty(n = 1)) +
    xlab("Time") +
    ylab(expression(S[y](t))) +
    labs(colour = "Set") +
    theme_light()
    
  p_IS <- data %>%
    ggplot(aes(x = block / n_blocks, y = IS, colour = set)) +
    scale_colour_discrete_sequential("ag_GrnYl") +
    geom_line() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_pretty(n = 1)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    xlab("Time") +
    ylab(expression(IS[y](t))) +
    labs(colour = "Set") +
    theme_light()
    
  p_z <- data %>%
    ggplot(aes(x = block / n_blocks, y = z, colour = set)) +
    scale_colour_discrete_sequential("ag_GrnYl") +
    geom_line() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_pretty(n = 1)) +
    xlab("Time") +
    ylab(expression(Z[y](t))) +
    labs(colour = "Set") +
    theme_light()
    
  p_cv <- data %>%
    pivot_longer(cols = c(ks, cm), names_to = "type", values_to = "statistic") %>%
    mutate(cv = case_when(type == "cm" ~ 0.1939, .default = 0.8135)) %>%
    ggplot(aes(x = block / n_blocks, y = statistic, colour = set)) +
    facet_grid(type ~ ., scales = "free", labeller = labeller(.default = toupper)) +
    scale_colour_discrete_sequential("ag_GrnYl") +
    geom_line() +
    geom_hline(aes(yintercept = cv), colour = "black", linetype = "dashed") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = breaks_pretty(n = 1)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1), add = c(0, 0)), breaks = breaks_pretty(n = 3)) +
    xlab("Time") +
    ylab(expression("Norm of" ~ Z[y]*(t))) +
    labs(colour = "Set") +
    theme_light()

  
  p <- ggarrange(p_IS, p_z, p_cv, nrow = format[1], ncol = format[2], common.legend = TRUE, legend = "none")
  return(p)
  
}

