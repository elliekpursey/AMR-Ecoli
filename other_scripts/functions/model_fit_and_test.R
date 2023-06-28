# function to call model, test fit, save estimates to csv, and save prediction
# plot to object

model_fit_and_test <- function(model_df, response, predictor, x_label, y_label){
  
  # create model formula
  
  form <- paste0(deparse(substitute(response)), ' ~ ', deparse(substitute(predictor)))
  
  predictor_string <- as.character(deparse(substitute(predictor)))
  
  # run model
  eval(substitute(model <- glmmTMB(formula(form),
                                   data=model_df, na.action = "na.fail", family = "binomial")))
  
  
  # define model name for output R object and write
  
  model_name <- paste0("model_objects/", deparse(substitute(model_df)),
                 "_",deparse(substitute(response)))
  
  saveRDS(model, model_name)
  
  # get model summary
  
  model_summary <- summary(model)
  
  print(model_summary)
  
  # test fit of model (output plots and statistics)
  
  print(testDispersion(model))
  simulationOutput <- simulateResiduals(fittedModel = model)
  plot(simulationOutput)
  
  plotResiduals(simulationOutput, model_df[[predictor_string]], quantreg = T)
  
  # extract model coefficients
  
  model_coeff <- model_summary$coefficients$cond
  coeff_filename <- paste0("subsampled_model_estimates/", deparse(substitute(model_df)),
                           "_",deparse(substitute(response)),".csv")
  write.csv(model_coeff, file = coeff_filename)
  
  # create prediction dataframe and plot, save to global env
  
  pred_df <- ggpredict(model, terms = c(predictor_string), type='fixed', ci.lvl = 0.95) 
  
  pred_plot <- ggplot(pred_df, aes(x = x, y = predicted))+ 
    theme_light() +
    theme(text = element_text(size=15), 
          legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_point(size=3, position = position_dodge(width=0.8)) +
    geom_errorbar(aes(x = x, ymin = conf.low, ymax = conf.high), size=0.5, position = position_dodge(width=0.8)) +  
    labs(x = x_label, y = y_label)
  
  plot_name <- paste0(deparse(substitute(model_df)),
                      "_",deparse(substitute(response)))
  print(plot_name)
  
  assign(plot_name, pred_plot, envir = .GlobalEnv)
  
}