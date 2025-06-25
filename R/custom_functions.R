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

# show conditions
show_conditions <- function(intervention_info){
  
  # Columns of interest
  target_cols <- c(
    "intervention_description",
    "intervention_selection_description",
    "control_selection_description"
  )
  
  intervention_info |> 
    select(any_of(target_cols)) |> 
    kbl()
  
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
    filter(effect == "fixed") %>% 
    mutate(
      # make sdt outcomes
      SDT_term = case_when(
        term == "(Intercept)" ~ "average c",
        term == "veracity_numeric" ~ "average d'",
        term == "condition_numeric" ~ "delta c",
        term == "veracity_numeric:condition_numeric" ~ "delta d'",
        TRUE ~ "Other"
      ), 
      # reverse c and delta c estimates and confidence intervals
      SDT_estimate = ifelse(term == "(Intercept)" | term == "condition_numeric", 
                            -1*estimate, estimate),
      SDT_conf.low = ifelse(term == "(Intercept)" | term == "condition_numeric", 
                            -1*conf.low, conf.low),
      SDT_conf.high = ifelse(term == "(Intercept)" | term == "condition_numeric", 
                            -1*conf.high, conf.high),
      sampling_variance = std.error^2
    ) 
  
  return(model)
  
}

# loop over experiments
run_loop <- function(data, filename){
  
  # only execute the following if the file does NOT exist
  if (!file.exists(filename)) {
    
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
    
    write_csv(results, filename)
    
    print(paste("Elapsed time: ", round(time[3]/60, digits = 2), " minutes"))
  }
  
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



