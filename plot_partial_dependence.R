gridMe <- function(x, length.out = 100){
  if(is.numeric(x)){
    seq(min(x),max(x),length.out = length.out)
  }
  else{
    
    rep(unique(x),length.out = length.out)
    
  }
  
}

smooth.spline2 <- function(formula, data, ...) { 
  mat <- model.frame(formula, data) 
  smooth.spline(mat[, 2], mat[, 1]) 
} 

predictdf.smooth.spline <- function(model, xseq, se, level) { 
  pred <- predict(model, xseq) 
  data.frame(x = xseq, y = pred$y) 
} 

plot_partial_dependence.Internal <- function(
  data,
  rf_model,
  r,
  n_bs_samples = 100,
  cov_set,
  n_grid_points
){
  
  print(head(data))
  data_list <- list()
  
  for(boot_sample in 1:n_bs_samples){
    
    data_list[[boot_sample]] <- list()
    
    #Resample data:
    boot_id <- sample(1:nrow(data),size = nrow(data), replace = TRUE)
    boot_data <- data[boot_id,]
    
    #Fit the model:
    formula <- as.formula(paste0(r," ~ ."))
    rf <- randomForest(
      formula = formula,
      data = boot_data,
      importance=FALSE,
      ntree=rf_model$ntree,
      mtry = rf_model$mtry)
    
    response_column <- which(colnames(data) == r)
    
    for(var_index in 1:length(cov_set)){
      #Check if numeric:
      IS_NUMERIC <- is.numeric(data[[cov_set[var_index]]])
      
      #Get predicted values:
      grid_data <- boot_data[rep(1:nrow(boot_data),times = n_grid_points),-response_column]
      
      #Grid:
      if(IS_NUMERIC){
        #Grid on original scale then convert to transformed scale:
        orig_data_grid <- gridMe(data[[cov_set[var_index]]],length.out = n_grid_points)
        grid_data[,cov_set[var_index]] <- rep(
          orig_data_grid,
          each = nrow(boot_data))
      }else{
        grid_data[,cov_set[[var_index]]] <- rep(
          gridMe(data[[cov_set[var_index]]],length.out = n_grid_points),
          each = nrow(boot_data))
      }
      
      #Predict:
      grid_data$prediction <- predict(rf,newdata = grid_data)
      grid_data <- grid_data[,c(cov_set[[var_index]],"prediction")]
      #Summarise:
      data_list[[boot_sample]][[var_index]] <- grid_data %>% group_by(UQ(as.name(cov_set[[var_index]]))) %>%
        summarise(med = mean(prediction)) %>% ungroup()
      
    }
  }
  return(data_list)
}

# source("functions/gridMe.R")

plot_partial_dependence <- function(
  data,
  rf_model,
  r,
  R = 100,
  n_cores = 5,
  cov_set,
  n_grid_points,
  inverse_covariate_transforms = NULL
){
  
  temp_func <- function(x){
    plot_partial_dependence.Internal(
      data = data,
      rf_model = rf_model,
      r = r,
      n_bs_samples = ceiling(R/n_cores),
      cov_set = cov_set,
      n_grid_points = n_grid_points
    )
  }
  
  t1 <- proc.time()
  cl<-makeCluster(n_cores)
  clusterExport(cl,
                varlist = c(
                  "plot_partial_dependence.Internal","R","n_cores",
                  "data","cov_set","n_grid_points","gridMe","r","rf_model"),envir = environment())
  clusterCall(cl, function(){library(randomForest);library(tidyverse)})
  bs_samples <- clusterApply(cl, 1:n_cores,temp_func)
  stopCluster(cl)
  t2 <- proc.time() - t1
  
  bs_samples <- do.call(args = bs_samples,what = "c")
  
  #Reorder into list of bs samples for each variable.
  var_bs_results_list <- list()
  for(var in 1:length(cov_set)){
    temp <- lapply(bs_samples,function(x) x[[var]])
    temp <- do.call("rbind", temp)
    
    temp$bs_sample <- rep(1:(ceiling(R/n_cores) * n_cores),each = length(unique(temp[,1][[1]])))
    var_bs_results_list[[var]] <- temp
  }
  names(var_bs_results_list) <- cov_set
  plot_list <- list()
  for(k in 1:length(cov_set)){
    IS_NUMERIC <- is.numeric(var_bs_results_list[[k]][,1][[1]])
    if(IS_NUMERIC){
      ci_data <- 
        var_bs_results_list[[k]] %>% group_by(bs_sample) %>% nest() %>%
        mutate(
          model = map(
            .x = data, 
            .f = function(x) smooth.spline2(data = x[,1:2], formula = med ~ .))) %>%
        mutate(
          xseq = map(
            .x = data, 
            .f = function(x) gridMe(x[[1]],50))) %>%
        mutate(predictions = map2(
          .x = model,
          .y = xseq,
          .f = function(x,y) predictdf.smooth.spline(model = x, xseq = y,se = FALSE, level = NULL)
        )) %>% 
        dplyr::select(bs_sample,predictions) %>% unnest(cols = c(predictions)) %>% 
        group_by(x) %>%
        summarise(med = median(y), lower = quantile(x = y, prob = 0.025), upper = quantile(x = y, prob = 0.975))
      
      # if(!is.null(inverse_covariate_transforms){
      #   ci_data$x <- (inverse_covariate_transforms[[cov_set[[k]]]])(ci_data$x)
      # }
      
      plot_list[[k]] <- 
        ggplot(data = ci_data) + 
        geom_line(aes(x = x, y = med)) + 
        geom_ribbon(aes(x = x,ymin = lower, ymax = upper),alpha = 0.4) +
        labs(y = NULL, x = cov_set[k]) +
        geom_rug(data = data,aes_string(x = cov_set[k]))
    }else{
      
      ci_data <- 
        var_bs_results_list[[k]] %>% group_by_(.dots = c(names(var_bs_results_list)[[k]])) %>% 
        summarise(
          mean = mean(med),
          medi = median(med),
          lower = quantile(
            x = med, prob = 0.025), 
          upper = quantile(x = med, prob = 0.975)
        )
      
      plot_list[[k]] <- ggplot(data = ci_data) + 
        geom_point(aes_string(x = names(var_bs_results_list)[k], y = "medi")) + 
        geom_errorbar(aes_string(x = names(var_bs_results_list)[k],ymin = "lower", ymax = "upper"),alpha = 0.4) +
        labs(y = NULL, x = cov_set[k])
      
    }
  }
  return(plot_list)
}