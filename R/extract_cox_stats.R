
# Function to extract Cox proportional hazards statistics for a given variable, endpoint and threshold quantile
extract_cox_stats <- function(data, endpoint, variable, selected_quantile) {
  time_col <- paste0(endpoint, "_years")
  
  # remove patients with missing endpoint information
  data = data %>%
    filter(!is.na(.data[[time_col]]),
           !is.na(.data[[time_col]]))
  
  # estimate best cutpoint
  if(selected_quantile == "best") {
    cutpoint <- surv_cutpoint(data = data,
                              time = time_col, event = endpoint,
                              variables = variable,
                              minprop = 0.1)[["cutpoint"]][["cutpoint"]]
    # estimate the empirical CDF and find corresponding quantile of the best cutpoint
    selected_quantile <- ecdf(data[[variable]])(cutpoint)
    
    # validate provided quantile to be between 0 and 1
  } else if (is.numeric(selected_quantile) & selected_quantile > 0 & selected_quantile < 1) {
    # determine threhsold based on provided quantile
    cutpoint <- quantile(data[[variable]], selected_quantile)
    
  } else {
    stop("selected_quantile must be 'best' or a numeric value between 0 and 1")
  }
  
  cox_formula <- glue("Surv({endpoint}_years, {endpoint}) ~ {variable} > {cutpoint}") %>% as.formula()
  cox_res <- coxph(cox_formula, data = data) %>%
    tidy(conf.int = TRUE, exponentiate = T) %>%
    # format output
    transmute(variable = variable,
              endpoint = endpoint,
              threshold = cutpoint,
              quantile = selected_quantile,
              HR = estimate,
              CI = glue("{round(conf.low, 2)} - {round(conf.high, 2)}"),
              p.value = p.value
    ) 
  
  return(cox_res)
}
