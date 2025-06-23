# Sensitivity simulations

# For detailed explanations, see `simulate_data.Rmd`

library("lme4")        # model specification / estimation
library("lmerTest")    # provides p-values in the output
library("tidyverse")   # data wrangling and visualisation
library("afex")        # anova and deriving p-values from lmer
library("broom")       # extracting data from model fits 
library("broom.mixed") # using tidy() for mixed models
library("metafor")     # doing mata analysis
library("kableExtra")  # for tables

# make sure were in the project directory
library(here)

# simulate a single sample
draw_sample <- function(
    # fixed effects
  beta_0,
  beta_v,
  beta_c,
  beta_cv,
  # random effects
  subj_0,
  subj_v,
  subj_rho_0v,
  # simulation related
  n_subj,
  n_fake,
  n_true,
  sigma,
  ...  # This allows additional parameters to be ignored
) { 
  
  # simulate a sample of items
  n_items <- n_fake + n_true
  
  items <- data.frame(
    item_id = seq_len(n_items),
    veracity = rep(c("fake", "true"), c(n_fake, n_true)),
    # get a numeric version of veracity that is effect-coded (i.e. not 0 vs. 1, 
    # but -0.5 and 0.5)
    veracity_effect_code = rep(c(-0.5, 0.5), c(n_fake, n_true))
  )
  
  # simulate a sample of subjects
  
  # calculate random intercept / random slope covariance
  covar <- subj_rho_0v * subj_0 * subj_v
  
  # put values into variance-covariance matrix
  cov_mx  <- matrix(
    c(subj_0^2, covar,
      covar,   subj_v^2),
    nrow = 2, byrow = TRUE)
  
  # generate the by-subject random effects
  subject_rfx <- MASS::mvrnorm(n = n_subj,
                               mu = c(SU_0 = 0, SU_v = 0),
                               Sigma = cov_mx)
  
  # combine with subject IDs
  subjects <- data.frame(subject_id = seq_len(n_subj),
                         subject_rfx)
  
  # cross subject and item IDs and calculate accuracy
  crossing(subjects, items)  %>%
    mutate(e_s = rnorm(nrow(.), mean = 0, sd = sigma))
}

# simulate an experiment
draw_experiment <- function(n_max_conditions = NULL, ...) {
  
  # generate random number of conditions
  possible_n_conditions <- seq(from = 2, to = n_max_conditions, by = 1)
  n_conditions <- sample(possible_n_conditions, size = 1)
  
  n_interventions <- n_conditions - 1 # we assume always one control group
  
  # draw control condition
  control <- draw_sample(...) %>% 
    mutate(condition = "control")
  
  # draw interventions
  interventions <-  1:n_interventions %>%
    map_dfr(function(x) {
      
      # To keep track of progress
      print(paste("drawing intervention number ", x))
      
      single_intervention <- draw_sample(...) %>% 
        mutate(condition = "intervention", 
               intervention_id = x)
      
      return(single_intervention)
    })
  
  # combine control and interventions data
  experiment <- bind_rows(control, interventions)
  
} 

# simulate multiple experiments
draw_multiple_experiments <- function(
    n_max_experiments, 
    exp_0, # by-experiment intercept sd
    exp_v, # by-experiment sensitivity (d prime) sd
    exp_c, # by-experiment - delta response bias (i.e. - delta c) sd
    exp_cv, # by-experiment delta sensitivity (i.e. delta d prime) sd
    exp_rho_0v, # correlation between intercept (- average response bias) and sensitivity (d prime)
    exp_rho_0c, # correlation between intercept and - delta response bias (i.e. - delta c)
    exp_rho_0cv, # correlation between intercept (- average response bias) and delta sensitivity (i.e. delta d prime)
    exp_rho_vc, # correlation between sensitivity (d prime) and - delta response bias (i.e. - delta c)
    exp_rho_vcv, # correlation between sensitivity (d prime) and delta sensitivity (i.e. delta d prime)
    exp_rho_ccv, # correlation between - delta response bias (i.e. - delta c) and delta sensitivity (i.e. delta d prime)
    ...) {
  
  # generate random number of experiments
  possible_n_experiments <- seq(from = 1, to = n_max_experiments, by = 1)
  n_experiments <- sample(possible_n_experiments, size = 1)
  
  # draw experiments
  experiments <-  1:n_experiments %>%
    map_dfr(function(x) {
      
      # To keep track of progress
      print(paste("drawing experiment number ", x))
      
      single_experiment <- draw_experiment(...) %>% 
        mutate(experiment_id = x)
      
      return(single_experiment)
    })
  
  # Calculate covariances
  covar_0v_exp <- exp_rho_0v * exp_0 * exp_v
  covar_0c_exp <- exp_rho_0c * exp_0 * exp_c
  covar_0cv_exp <- exp_rho_0cv * exp_0 * exp_cv
  covar_vc_exp <- exp_rho_vc * exp_v * exp_c
  covar_vcv_exp <- exp_rho_vcv * exp_v * exp_cv
  covar_cc_exp <- exp_rho_ccv * exp_c * exp_cv
  
  # Create the variance-covariance matrix
  cov_mx_exp <- matrix(
    c(exp_0^2, covar_0v_exp, covar_0c_exp, covar_0cv_exp,
      covar_0v_exp, exp_v^2, covar_vc_exp, covar_vcv_exp,
      covar_0c_exp, covar_vc_exp, exp_c^2, covar_cc_exp,
      covar_0cv_exp, covar_vcv_exp, covar_cc_exp, exp_cv^2),
    nrow = 4, byrow = TRUE
  )
  
  # Generate the by-experiment random effects
  exp_rfx <- MASS::mvrnorm(n = n_experiments,
                           mu = c(EXP_0 = 0, EXP_v = 0, EXP_c = 0, EXP_cv = 0),
                           Sigma = cov_mx_exp)
  
  # Check the value of n_experiments
  if (n_experiments == 1) {
    # if n == 1, mvnorm seems to return a vector, which data.frame then turns 
    # into long format data (instead of wide-format when n_experiments > 1)
    # Reshape to wide format
    exp_rfx <- data.frame(exp_rfx) %>%
      mutate(effect = rownames(.)) %>% 
      pivot_wider(names_from = effect, values_from = "exp_rfx")
  }
  
  # Combine with experiment IDs
  exp_rfx <- data.frame(experiment_id = seq_len(n_experiments), exp_rfx)
  
  # add random effects to experiment data
  experiments <- left_join(experiments, exp_rfx, by = "experiment_id")
  
  return(experiments)
}

# simulate papers
draw_papers <- function(
    paper_0, # analogous to experiment random effects
    paper_v,
    paper_c,
    paper_cv,
    paper_rho_0v,
    paper_rho_0c,
    paper_rho_0cv,
    paper_rho_vc,
    paper_rho_vcv, 
    paper_rho_ccv,
    n_papers = 10, # number of papers in the meta analysis
    ...) {
  
  # draw papers
  papers <-  1:n_papers %>%
    map_dfr(function(x) {
      
      # To keep track of progress
      print(paste("drawing paper number ", x))
      
      single_paper <- do.call(draw_multiple_experiments, parameters) %>% 
        mutate(paper_id = x)
      
      return(single_paper)
    })
  
  # Calculate covariances
  covar_0v_paper <- paper_rho_0v * paper_0 * paper_v
  covar_0c_paper <- paper_rho_0c * paper_0 * paper_c
  covar_0cv_paper <- paper_rho_0cv * paper_0 * paper_cv
  covar_vc_paper <- paper_rho_vc * paper_v * paper_c
  covar_vcv_paper <- paper_rho_vcv * paper_v * paper_cv
  covar_cc_paper <- paper_rho_ccv * paper_c * paper_cv
  
  # Create the variance-covariance matrix
  cov_mx_paper <- matrix(
    c(paper_0^2, covar_0v_paper, covar_0c_paper, covar_0cv_paper,
      covar_0v_paper, paper_v^2, covar_vc_paper, covar_vcv_paper,
      covar_0c_paper, covar_vc_paper, paper_c^2, covar_cc_paper,
      covar_0cv_paper, covar_vcv_paper, covar_cc_paper, paper_cv^2),
    nrow = 4, byrow = TRUE
  )
  
  # Generate the by-paper random effects
  paper_rfx <- MASS::mvrnorm(n = n_papers,
                             mu = c(PAP_0 = 0, PAP_v = 0, PAP_c = 0, PAP_cv = 0),
                             Sigma = cov_mx_paper)
  
  # Combine with paper IDs
  paper_rfx <- data.frame(paper_id = seq_len(n_papers), paper_rfx)
  
  # add random effects to paper data
  papers <- left_join(papers, paper_rfx) %>% 
    mutate(    
      # unique experiment identifier
      unique_experiment_id = paste(paper_id, experiment_id, sep = "_"), 
      # unique participant identifier
      unique_subject_id = paste0(paper_id, "_", experiment_id, intervention_id, "_", subject_id), 
      # unique intervention identifier
      unique_intervention_id = paste(paper_id, experiment_id, intervention_id, sep = "_"))
}

# calculate accuracy outcomes
calculate_accuracy <- function(
    data, 
    # fixed effects
    beta_0, # intercept (- average response bias )
    beta_v, # average sensitivity (d prime) 
    beta_c, # - delta response bias (i.e. - delta c)
    beta_cv, # delta sensitivity (i.e. delta d prime)
    ...
){
  
  data <- data %>% 
    mutate(
      # make numeric helper variables using deviation coding
      veracity_numeric = ifelse(veracity == "true", 0.5, -0.5),
      condition_numeric = ifelse(condition == "intervention",  0.5, -0.5),
      # calculate z-value accuracy
      accuracy_z_value = 
        (beta_0 + SU_0 + EXP_0 + PAP_0) +
        (beta_v + SU_v + EXP_v + PAP_v)*veracity_numeric +
        (beta_c + EXP_c + PAP_c )*condition_numeric+
        (beta_cv + EXP_cv + PAP_cv)*(condition_numeric*veracity_numeric) +
        e_s,
      # use the probit function (i.e. cumulative distribution function of the standard normal distribution)
      # to optain probabilities from the z-values
      accuracy_probability = pnorm(accuracy_z_value), 
      # generate a binary accuracy response by sampling from a bernoulli distribution
      # (i.e. binomial distribution with one trial)
      accuracy = rbinom(nrow(.), 1, accuracy_probability)
    ) %>% 
    # remove random effect variables
    select(-matches("SU|EXP|PAP", ignore.case = FALSE))
}

# final data simulation function
simulate_data <- function(...) {
  papers <- draw_papers(...)
  
  final_data <- calculate_accuracy(papers, ...)
  
  return(final_data)
}

# calculate model function
calculate_SDT_model <- function(data) {
  
  time <- system.time({
    model <- glmer(accuracy ~ veracity_numeric + condition_numeric + 
                     veracity_numeric*condition_numeric +
                     (1 + veracity_numeric | unique_subject_id),
                   data = data, 
                   family = binomial(link = "probit")
    )
  })
  
  time <- round(time[3]/60, digits = 2)
  
  # get a tidy version
  model <- tidy(model, conf.int = TRUE) %>% 
    # add time
    mutate(time_minutes = time)
  
  # give nicer names in SDT terminology to estimates (! and reverse estimates for response bias !)
  model <- model %>% 
    mutate(
      # reverse c and delta c estimates
      estimate = ifelse(term == "(Intercept)" | term == "condition_numeric", 
                        -1*estimate, estimate),
      term = case_when(term == "(Intercept)" ~ "average response bias (c)", 
                       term == "veracity_numeric" ~ "average sensitivity (d')", 
                       term == "condition_numeric" ~ "delta c",
                       term == "veracity_numeric:condition_numeric" ~ "delta d'",
                       TRUE ~ term
      ),
      sampling_variance = std.error^2
    ) 
  
  return(model)
  
}

# loop over experiments
run_SDT_model_loop <- function(data){
  
  # make a vector with all unique experiment ids
  experiments <- data %>% 
    distinct(unique_experiment_id) %>% 
    # slice(1:3) %>% # to test loop
    pull()
  
  time <- system.time({
    
    # run one model per experiment and store the results in a common data frame
    results <- experiments %>%
      map_dfr(function(x) {
        
        # restrict data to only the respective experiment
        experiment <- data %>% filter(unique_experiment_id == x)
        
        # extract paper id
        paper_id <- unique(experiment$paper_id)
        
        # To keep track of progress
        print(paste("calculating model for experiment ", x))
        
        model_experiment <- calculate_SDT_model(experiment) %>%
          mutate(unique_experiment_id = x,
                 paper_id = paper_id)
        
        return(model_experiment)
      })
  })
  
  print(paste("Elapsed time: ", round(time[3]/60, digits = 2), " minutes"))
  
  return(results)
}

# Function to calculate meta models
calculate_meta_model <- function(data, yi, vi, robust = TRUE) {
  
  # provide metafor compatible names
  metafor_data <- data %>% 
    rename(yi = {{yi}}, 
           vi = {{vi}})
  
  # Multilevel random effect model for accuracy
  model <-  metafor::rma.mv(yi, vi, random = ~ 1 | paper_id / unique_experiment_id, 
                            data = metafor_data)
  
  return(model)
  
  if(robust == TRUE) {
    # with robust standard errors clustered at the paper level 
    robust_model <- robust(model, cluster = data$paper_id)
    
    return(robust_model)
  }
}

# final function for running all models 
get_meta_estimates <- function(data){
  
  # generate outcome data for all experiments
  outcome_data <- run_SDT_model_loop(data)
  
  # calculate meta model for delta dprime
  delta_dprime <- calculate_meta_model(data = outcome_data  %>% 
                                         filter(term == "delta d'"), 
                                       yi = estimate, 
                                       vi = sampling_variance, 
                                       robust = TRUE) %>% 
    tidy() %>% 
    mutate(term = ifelse(term == "overall", "delta d'", NA))
  
  # calculate meta model for delta c
  delta_c <- calculate_meta_model(data = outcome_data  %>% 
                                    filter(term == "delta c"), 
                                  yi = estimate, 
                                  vi = sampling_variance, 
                                  robust = TRUE) %>% 
    tidy() %>% 
    mutate(term = ifelse(term == "overall", "delta c", NA))
  
  estimates <- bind_rows(delta_dprime, delta_c)
  
  return(estimates)
  
}

# repeat the process and store some information on it (e.g. time taken, number of iteration)
iterate <- function(iterations, ...) {
  
  # create data frame with model results for generated samples
  result <- 1:iterations %>%
    purrr::map_df(function(x) {
      
      # Measure time for each iteration
      iteration_time <- system.time({
        # generate data
        participant_level_data <- simulate_data(...)
        
        # run models on the data
        estimates <- get_meta_estimates(participant_level_data)
      })
      
      # Add iteration number and time in minutes to estimates
      estimates <- estimates %>%
        mutate(iteration = x,
               time_taken_minutes = iteration_time[3] / 60)  # Convert time to minutes
      
      # To keep track of progress
      if (x %% 2 == 0) {print(paste("iteration number ", x))}
      
      return(estimates)
    })
  
  return(result)
}

# run sensitivity analysis
sensitivity_analysis <- function(file_name, effect_sizes_to_try, iterations, ...) {
  
  # make full file name
  full_file_name <- paste0("data/simulations/", file_name)
  
  # only run analysis if a file with that name does not yet exists
  if (!file.exists(full_file_name)) {
    
    sensitivity_data <- pmap_df(effect_sizes_to_try, function(delta_d_prime, delta_c) {
      
      # To keep track of progress
      print(paste("test for delta c = ", delta_c, "; delta d' = ", delta_d_prime))
      
      # calculate the results for a specific combination of parameters
      combination_result <-  iterate(iterations = iterations, 
                                     beta_cv = delta_d_prime, 
                                     # make sure to reverse the data generating parameter for response bias
                                     beta_c = -delta_c, 
                                     ...
      )
      
      # add combination of parameters to data frame
      combination_result <- combination_result %>% 
        mutate(parameter_delta_d_prime = delta_d_prime, 
               parameter_delta_c = delta_c)
      
    })
    
    write_csv(sensitivity_data, full_file_name)
    
    return(sensitivity_data)
  }
}

# Set parameters
parameters <- list(
  # fixed effects
  beta_0   = -0.5, # intercept (- average response bias )
  beta_v   = 1, # average sensitivity (d prime) 
  beta_c   = - 0.1, # - delta response bias (i.e. - delta c)
  beta_cv   = 0.1, # delta sensitivity (i.e. delta d prime)
  # random effects
  # subjects
  subj_0   = 0.1, # by-subject intercept sd
  subj_v   = 0.1, # by-subject sensitivity (d prime) sd
  subj_rho_0v = 0, # correlation between intercept (- average response bias) and sensitivity (d prime)
  # experiments
  exp_0   = 0.1, # by-experiment intercept sd
  exp_v   = 0.1, # by-experiment sensitivity (d prime) sd
  exp_c   = 0.1, # by-experiment - delta response bias (i.e. - delta c) sd
  exp_cv   = 0.1, # by-experiment delta sensitivity (i.e. delta d prime) sd
  exp_rho_0v = 0, # correlation between intercept (- average response bias) and sensitivity (d prime)
  exp_rho_0c = 0, # correlation between intercept and - delta response bias (i.e. - delta c)
  exp_rho_0cv = 0, # correlation between intercept (- average response bias) and delta sensitivity (i.e. delta d prime)
  exp_rho_vc = 0, # correlation between sensitivity (d prime) and - delta response bias (i.e. - delta c)
  exp_rho_vcv = 0, # correlation between sensitivity (d prime) and delta sensitivity (i.e. delta d prime)
  exp_rho_ccv = 0, # correlation between - delta response bias (i.e. - delta c) and delta sensitivity (i.e. delta d prime)
  # papers
  paper_0   = 0.1, # analogous to experiment random effects
  paper_v   = 0.1, 
  paper_c   = 0.1, 
  paper_cv   = 0.1, 
  paper_rho_0v = 0, 
  paper_rho_0c = 0, 
  paper_rho_0cv = 0, 
  paper_rho_vc = 0, 
  paper_rho_vcv = 0, 
  paper_rho_ccv = 0, 
  sigma    = 0.5, # residual (error) sd (i.e. all variation that the model cannot account for)
  # simulation related
  n_subj  = 100, # number of subjects for each experimental arm (i.e. an experiment with a control and one intervention group has n = 200)
  n_fake  =  5, # number of fake news items 
  n_true  =  5, # number of true news items
  n_max_conditions = 4, # maximum number of possible conditions for a single experiment
  n_max_experiments = 4, # maximum number of possible experiments within a single paper
  n_papers = 10 # number of papers in the meta analysis
)

# Remove initial parameters from list since we want to replace it
# with the respective effect size combination every time.
parameters[["beta_cv"]] <- NULL
parameters[["beta_c"]] <- NULL

effect_sizes_to_try <- crossing(delta_d_prime = c(0.2, 0.5, 0.8),
                                delta_c = c(0.2, 0.5, 0.8))

do.call(sensitivity_analysis, c(parameters,
                                list(file_name = "sensitivity_analysis.csv",
                                     effect_sizes_to_try = effect_sizes_to_try,
                                     iterations = 100)
)
)

sensitivity_data <- read_csv("data/simulations/sensitivity_analysis.csv")












