# round numbers from models
rounded_numbers <- function(x) mutate_if(x, is.numeric, round, 3)

# for cleaning data from papers
save_data <- function(data) {
  filename_base <- deparse(substitute(data))
  
  # Save as CSV
  readr::write_csv(data, paste0("cleaned_", filename_base, ".csv"))
  
  # Save as RDS
  saveRDS(data, paste0("cleaned_", filename_base, ".rds"))
}


# get model outputs ready for inline text reporting
text_ready <- function(model_output) {
  
  result <- tidy(model_output, conf.int = TRUE) %>% 
    filter(effect == "fixed") %>% 
    # report p.value according to apa standards
    mutate(p.value = case_when(p.value < 0.001 ~ "< .001",
                               TRUE ~ sprintf("%.3f", p.value)
    )
    ) %>% 
    # all other terms
    rounded_numbers() %>% 
    mutate(
      term = ifelse(term == "(Intercept)", "intercept", term),
      ci = glue::glue("[{conf.low}, {conf.high}]"), 
      # if there is an interaction (recognized by ":"), name it just interaction
      term = str_replace(term, ".*:.*", "interaction")
    ) %>% 
    select(term, estimate, std.error, ci, p.value) %>% 
    split(.$term)
  
  return(result)
}

# clean data frame results for paired t-test
clean_t_test <- function(t.test) {
  
  t.test %>% 
    tidy() %>% 
    select(-c(method, alternative)) %>% 
    # report p.value according to apa standards
    mutate(p.value = case_when(p.value < 0.001 ~ "< .001",
                               TRUE ~ sprintf("p = %.3f", p.value)
    )
    ) %>% 
    # all other terms
    rounded_numbers() %>%
    mutate_all(~ {attributes(.) <- NULL; .})
  
}

# Plot results from simulation

# Plot competence distributions
plot_competence_distributions <- function(data) {
  # Your data frame
  d <- data %>% 
    group_by(competence) %>% 
    reframe(alpha = unique(alpha), beta = unique(beta))
  
  # Function to plot the different competence distributions 
  plot_competence_distribution <- function(alpha, beta, competence) {
    ggplot() +
      geom_function(
        fun = ~dbeta(.x, shape1 = alpha, shape2 = beta)
      ) +
      labs(title = paste0(competence,"\n", "alpha = ", alpha, ", beta = ", beta), 
           y = "Density", x = "Competence") +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      plot_theme +
      rremove("ylab") +  # Remove y-axis label
      rremove("xlab")   # Remove x-axis label
  }
  
  # Specify the values of alpha and beta from the data
  alpha_values <- d$alpha
  beta_values <- d$beta 
  competence_values <- d$competence
  
  # Create a list of ggplot objects
  plot_list <- pmap(list(alpha = alpha_values, beta = beta_values, competence = competence_values),
                    plot_competence_distribution)
  
  # Arrange and display the plots
  competence_plot <- ggarrange(plotlist = plot_list) %>% 
    annotate_figure(
                    left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 1)),
                    bottom = textGrob("Competence", gp = gpar(cex = 1)))
  
  
  return(competence_plot) 
}

# Function for splitting data along several variables (useful for inline reporting)
# taken from here: https://www.tjmahr.com/lists-knitr-secret-weapon/
super_split <- function(.data, ...) {
  dots <- rlang::enquos(...)
  for (var in seq_along(dots)) {
    var_name <- rlang::as_name(dots[[var]])
    .data <- purrr::map_depth(
      .x = .data,
      .depth = var - 1,
      .f = function(xs) split(xs, xs[var_name])
    )
  }
  .data
}


# Categorical scenario

# Vary competence
# group by constellations
plot_competence_vary_categorical <- function(data, variable = options, outcome = "Accuracy", 
                                 return = "simulations") {
  
  d <- data
  
  # Extract variable name as a string
  variable_name <- deparse(substitute(variable))
  
  # retrieve info from choice of variable
  if (variable_name == "options") {
    x_label = "Number of Choice Options"
    caption = paste0("Sample size of informants: ", d$sample[1],
                     "; Iterations per sample size: ", d$total_iterations[1], 
                     "; Population size per iteration: ", d$population[1])
  }
  if (variable_name == "sample") {
    x_label = "Size of informant groups"
    caption = paste0("Number of choice options: ", d$options[1],
                     "; Iterations per sample size: ", d$total_iterations[1], 
                     "; Population size per iteration: ", d$population[1])
  }
  
  # make constellation a factor with the right levels
  d$constellation <- fct_relevel(d$constellation, "minority(ies)", "dissensus", "majority", "consensus")
  
  d <- d %>% 
    group_by(constellation, {{variable}}, competence) %>% 
    summarise(competence_mean = mean(average_relative_competence),
              accuracy_mean = mean(average_accuracy),
              accuracy_sd = sd(average_accuracy),
              competence_sd = sd(average_relative_competence)
    )
  
  # Assign outcome variable
  if(outcome == "Accuracy") {
    d$outcome <- d$accuracy_mean
    d$outcome_sd <- d$accuracy_sd
  }
  
  if(outcome == "Competence") {
    d$outcome <- d$competence_mean
    d$outcome_sd <- d$competence_sd
  }
  
  # plot for accuracy
  plot_outcome <- ggplot(d,
                         aes(x = {{variable}}, y = outcome,
                             color = constellation)) +
    # geom_half_violin(position = position_nudge(x = -.2),
    #                  adjust=2, alpha = .5,
    #                  side = "l") +
    geom_pointrange(aes(ymin = outcome - outcome_sd, 
                        ymax = outcome + outcome_sd), 
                    size = 0.1) +
    geom_line() +
    # Add nice labels
    labs(x = x_label, y = outcome, caption = caption) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, 
                          limits = rev(levels(d$constellation)),
                          direction = -1
    ) +
    facet_wrap(~competence) +
    plot_theme + 
    theme(legend.position = "top")
  
  # Plot competence distributions
  # Your data frame
  d <- data %>% 
    group_by(competence) %>% 
    reframe(alpha = unique(alpha), beta = unique(beta))
  
  # Function to plot the different competence distributions 
  plot_competence_distributions <- function(alpha, beta) {
    ggplot() +
      geom_function(
        fun = ~dbeta(.x, shape1 = alpha, shape2 = beta)
      ) +
      labs(title = paste0("alpha = ", alpha, ", beta = ", beta), 
           y = "Density", x = "Competence") +
      plot_theme
  }
  
  # Specify the values of alpha and beta from the data
  alpha_values <- d$alpha
  beta_values <- d$beta 
  
  # Create a list of ggplot objects
  plot_list <- map2(alpha_values, beta_values, plot_competence_distributions)
  
  # Arrange and display the plots
  competence_plot <- ggarrange(plotlist = plot_list) %>%
    ggpubr::annotate_figure(top = ggpubr::text_grob("Competence distributions", 
                                                    face = "bold",
                                                    color = "darkgrey",
                                                    size = 14))
  
  if(return == "competence distributions") {
    return(competence_plot)
  } else {
    return(plot_outcome) 
  }
} 

# function used in participant vs. model comparisons
plot_results <- function(data, outcome = "everything") {
  
  d <- data
  
  # make constellation a factor with the right levels
  d$constellation <- fct_relevel(d$constellation, "minority(ies)", "dissensus", "majority", "consensus")
  
  # extract simulation info
  simulation_info <- d %>% summarize(across(c(population, sample, options, total_iterations), mean)) 
  
  caption = paste0("Number of choice options: ", simulation_info$options,
                   "; Sample size: ", simulation_info$sample,
                   "; Iterations : ", simulation_info$total_iterations, 
                   "; Population size per iteration: ", simulation_info$population)
  
  # make descriptive data
  # descriptive <- d %>% 
  #   group_by(constellation) %>% 
  #   summarise(count = sum(count)) %>% 
  #   mutate(rel_freq = round(count / sum(count), digits = 2)) %>% 
  #   gridExtra::tableGrob(rows = NULL)
  
  
  # plot for accuracy
  plot_accuracy <- ggplot(d,
                          aes(x = constellation, y = average_accuracy, fill = constellation)) +
    geom_half_violin(position = position_nudge(x = -.2),
                     adjust=2, alpha = .8,
                     side = "l") +
    stat_summary(fun = "mean", geom = "point", size = 1, shape = 21) +
    stat_summary(fun = "mean", geom = "line", size = 1, linetype = "dashed") +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
    # Add nice labels
    labs(x = "Convergence", y = "Accuracy") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1) +
    guides(fill = FALSE) +
    plot_theme + 
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    ylim(c(0, 1))
  
  # plot for competence
  plot_competence <- ggplot(d,
                            aes(x = constellation, y =  average_relative_competence, fill = constellation)) +
    geom_half_violin(position = position_nudge(x = -.2),
                     adjust=2, alpha = .8,
                     side = "l") +
    stat_summary(fun = "mean", geom = "point", size = 1, shape = 21) +
    stat_summary(fun = "mean", geom = "line", size = 1, linetype = "dashed") +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
    # Add nice labels
    labs(x = "Convergence", y = "Competence", caption = caption) +
    scale_fill_viridis_d(option = "plasma", begin = 0.1) +
    guides(fill = FALSE) +
    plot_theme + 
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    ylim(c(0, 1))
  
  # unite tables
  # tables <- gridExtra::grid.arrange(descriptive, simulation_info, ncol = 2)
  
  if (outcome == "everything") {
    
    return(plot_accuracy | plot_competence )
  }
  
  if (outcome == "accuracy") {
    
    return(plot_accuracy + labs(caption = caption))
  }
  
  if (outcome == "competence") {
    
    return(plot_competence)
  }
  
}

# main plot by size of relative majority
plot_competence_vary_relative_majority_categorical <- function(data, variable = options, outcome = "Accuracy") {
  
  d <- data
  
  # Extract variable name as a string
  variable_name <- deparse(substitute(variable))
  
  # retrieve info from choice of variable
  if (variable_name == "options") {
    caption = paste0("Sample size of informants: ", d$sample[1],
                     "; Iterations per sample size: ", d$total_iterations[1], 
                     "; Population size per iteration: ", d$population[1])
  }
  if (variable_name == "sample") {
    caption = paste0("Number of choice options: ", d$options[1],
                     "; Iterations per sample size: ", d$total_iterations[1], 
                     "; Population size per iteration: ", names(table(d$population))[1], "/", names(table(d$population))[2])
  }
  
  d <- d %>% 
    group_by(relative_majority, {{variable}}, competence) %>%
    summarise(competence_mean = mean(average_relative_competence),
              accuracy_mean = mean(average_accuracy),
              accuracy_sd = sd(average_accuracy),
              competence_sd = sd(average_relative_competence)
    )
  
  # Assign outcome variable
  if(outcome == "Accuracy") {
    d$outcome <- d$accuracy_mean
    d$outcome_sd <- d$accuracy_sd
  }
  
  if(outcome == "Competence") {
    d$outcome <- d$competence_mean
    d$outcome_sd <- d$competence_sd
  }
  
  plot_outcome <- ggplot(d, 
                         aes(x = relative_majority, y = outcome,
                             color = as.factor({{variable}}))) +
    geom_point() +
    geom_line() +
    # Add nice labels
    labs(x = "Convergence \n (share of votes for an option)", 
         y = outcome, caption = caption, 
         color = paste0("n ", variable_name)) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, 
                          #limits = rev(levels(d$constellation)),
                          direction = -1
    ) +
    facet_wrap(~competence) +
    plot_theme + 
    theme(legend.position = "top")
  
  return(plot_outcome) 
  
} 

# Numerical scenario

# plot a single one sample size
plot_competence_single_sample <- function(data, outcome = "accuracy") {
  
  d <- data
  
  # extract simulation info
  caption = paste0("Sample size of informants: ", d$sample[1],
                   "; Iterations per sample size: ", d$total_iterations[1], 
                   "; Population size per iteration: ", d$population[1])
  
  
  plot_accuracy <- ggplot(d, aes(x = convergence, y = accuracy_mean)) +
    geom_hex(alpha = 0.9) +
    #geom_smooth(method = "lm", se = FALSE) +
    stat_smooth(geom='line', method = "lm", alpha=0.7, se=FALSE, color = "darkgrey") +
    # Add nice labels
    labs(x = "Convergence \n(SDs of samples)", 
         y = "Accuracy \n(squared distance to true mean)",
         caption = caption) +
    scale_fill_viridis_c(option = "plasma", begin = 0.1) +
    facet_wrap(~competence) +
    plot_theme + 
    theme(legend.position = "top", 
          legend.key.width = unit(1.5, "cm"), 
          legend.key.height = unit(0.3, "cm"))
  
  plot_competence <- ggplot(d, aes(x = convergence, y = competence_mean)) +
    geom_hex(alpha = 0.9) +
    #geom_smooth(method = "lm", se = FALSE) +
    stat_smooth(geom='line', method = "lm", alpha=0.7, se=FALSE, color = "darkgrey") +
    # Add nice labels
    labs(x = "Convergence \n(SDs of samples)", 
         y = "Competence \n(SD for data generating function)",
         caption = caption) +
    scale_fill_viridis_c(option = "plasma", begin = 0.1) +
    facet_wrap(~competence) +
    plot_theme + 
    theme(legend.position = "top", 
          legend.key.width = unit(1.5, "cm"), 
          legend.key.height = unit(0.3, "cm"))
  
  if (outcome == "accuracy") {
    
    return(plot_accuracy)
  }
  
  if (outcome == "competence") {
    
    return(plot_competence)
  }
}

# plot a list by sample size
plot_competence_vary_numeric <- function(data, ...) {
  
  # get different samples as vector
  samples <- data %>% 
    reframe(unique(as.character(sample))) %>% pull()
  
  
  # empty list
  plot_list <- list()
  
  # make plots fro each moderator
  for (i in samples) {
    
    plot_data <- data %>%
      filter(sample == i)
    
    sample_plot <- plot_competence_single_sample(data = plot_data, ...)
    
    plot_list[[i]] <- sample_plot
  }
  
  return(plot_list)
}

# plot continuous simulations for cogsci paper
plot_continuous_cogsci <- function(data, outcome = "accuracy") {
  
  d <- data
  
  # extract simulation info
  caption = paste0("Iterations per sample size: ", d$total_iterations[1], 
                   "; Population size per iteration: ", d$population[1])
  
  
  if(outcome == "accuracy") {
    
    d$outcome <- d$accuracy_mean
    
    outcome_label <- "Accuracy"
  }
  
  if (outcome == "competence") {
    
    d$outcome <- d$competence_mean
    
    outcome_label <- "Competence"
  }  
  
  plot <- ggplot(d, aes(x = convergence, y = outcome, color = as.factor(sample))) +
    #geom_smooth(method = "lm", se = FALSE) +
    stat_smooth(geom='line', method = "lm", alpha=0.7, se=FALSE) +
    # Add nice labels
    labs(x = "Convergence \n(SDs of samples)", 
         y = outcome_label,
         caption = caption, 
         color = "n sample") +
    scale_color_viridis_d(option = "plasma", begin = 0.1, 
                          #limits = rev(levels(d$constellation)),
                          direction = -1
    ) +
    facet_wrap(~competence) +
    plot_theme + 
    theme(legend.position = "top", 
          legend.key.width = unit(1.5, "cm"), 
          legend.key.height = unit(0.3, "cm"))

  return(plot)
}

# combined plot for all simulation scenarios

plot_continuous_combined <- function(data, outcome = "accuracy") {
  
  d <- data %>% 
    mutate(sample = factor(sample, levels = sample, labels = paste0("N groups = ", sample)))
  
  # Maybe later:
  # make a relative variable for the fill aesthetics 
  # (the default mode for geom_hex is to count N, but this differs so much across simulations that legend is not meaningful)
  # rel_frequency_in_cell <- data %>% 
  #   group_by(competence, sample, accuracy, ) %>% 
  #   summarise(n = n(), 
  #             rel.)
  
  
  # extract simulation info
  caption = paste0("Iterations per sample size: ", d$total_iterations[1], 
                   "; Population size per iteration: ", d$population[1])
  
  
  plot_accuracy <- ggplot(d, aes(x = convergence, y = accuracy_mean)) +
    geom_hex(alpha = 0.9) +
    #geom_smooth(method = "lm", se = FALSE) +
    stat_smooth(geom='line', method = "lm", alpha=0.7, se=FALSE, color = "darkgrey") +
    # Add nice labels
    labs(x = "Convergence \n(SDs of samples)", 
         y = "Accuracy \n(squared distance to true mean)",
         caption = caption) +
    scale_fill_viridis_c(option = "plasma", begin = 0.1) +
    facet_grid(sample ~ competence) +
    plot_theme + 
    theme(legend.position = "top", 
          legend.key.width = unit(1.5, "cm"), 
          legend.key.height = unit(0.3, "cm"))
  
  plot_competence <- ggplot(d, aes(x = convergence, y = competence_mean)) +
    geom_hex(alpha = 0.9) +
    #geom_smooth(method = "lm", se = FALSE) +
    stat_smooth(geom='line', method = "lm", alpha=0.7, se=FALSE, color = "darkgrey") +
    # Add nice labels
    labs(x = "Convergence \n(SDs of samples)", 
         y = "Competence \n(SD for data generating function)",
         caption = caption) +
    scale_fill_viridis_c(option = "plasma", begin = 0.1) +
    facet_grid(sample ~ competence) +
    plot_theme + 
    theme(legend.position = "top", 
          legend.key.width = unit(1.5, "cm"), 
          legend.key.height = unit(0.3, "cm"))
  
  if (outcome == "accuracy") {
    
    return(plot_accuracy)
  }
  
  if (outcome == "competence") {
    
    return(plot_competence)
  }
}


# old functions that might still be useful at some point

# Plot competence distributions
OLD_plot_competence_distributions <- function(data) {
  # Your data frame
  d <- data %>% 
    group_by(competence) %>% 
    reframe(alpha = unique(alpha), beta = unique(beta))
  
  # Function to plot the different competence distributions 
  plot_competence_distribution <- function(alpha, beta) {
    ggplot() +
      geom_function(
        fun = ~dbeta(.x, shape1 = alpha, shape2 = beta)
      ) +
      labs(title = paste0("alpha = ", alpha, ", beta = ", beta), 
           y = "Density", x = "Competence") +
      plot_theme
  }
  
  # Specify the values of alpha and beta from the data
  alpha_values <- d$alpha
  beta_values <- d$beta 
  
  # Create a list of ggplot objects
  plot_list <- map2(alpha_values, beta_values, plot_competence_distribution)
  
  # Arrange and display the plots
  competence_plot <- ggarrange(plotlist = plot_list) %>%
    ggpubr::annotate_figure(top = ggpubr::text_grob("Competence distributions", 
                                                    face = "bold",
                                                    color = "darkgrey",
                                                    size = 14))
  
  
  return(competence_plot) 
}









