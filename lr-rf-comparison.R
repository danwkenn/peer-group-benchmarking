
{
lm_bootresids_OOB <- bootstrap_residuals(
                    orig_data = data_list[[r]],
                    caret_fit = caret_lm_fit,
                    n_resamples = 1000,
                    ncores = 5,
                    verbose = TRUE,
                    force_OOB = TRUE,
                    strata = clusters
                    )

lm_bootresids_IB <- bootstrap_residuals(
  orig_data = data_list[[r]],
  caret_fit = caret_lm_fit,
  n_resamples = 1000,
  ncores = 5,
  verbose = TRUE,
  force_OOB = FALSE,
  strata = clusters
)

plot(density(c(lm_bootresids_OOB),na.rm = TRUE))
lines(density(lm_bootresids_IB),col = "red")
# OOB matters for linear regression too.


plot(density(c(boot_resids1),na.rm = TRUE))
lines(density(lm_bootresids_OOB,na.rm = TRUE),col = "red")
sd(boot_resids1,na.rm = TRUE)
sd(lm_bootresids,na.rm = TRUE)

comparison_data <- data.frame(comparison=character(),RF = numeric(), LR = numeric())
temp <- 
  data.frame(
    comparison = "Mean",
    RF = apply(boot_resids1,1, mean,na.rm = TRUE),
    LR = apply(lm_bootresids_OOB,1, mean,na.rm = TRUE)
  )
summary(lm(data = temp, I(abs(LR)-1*abs(RF)) ~ (abs(RF))))
comparison_data <- rbind(
  comparison_data,
  temp
)
temp2 <- 
  data.frame(
    comparison = "SD",
    RF = apply(boot_resids1,1, sd,na.rm = TRUE),
    LR = apply(lm_bootresids_OOB,1, sd,na.rm = TRUE)
  )
comparison_data <- rbind(
  comparison_data,
  temp2
)
temp <- 
  data.frame(
    comparison = "median",
    RF = apply(boot_resids1,1, median,na.rm = TRUE),
    LR = apply(lm_bootresids_OOB,1, median,na.rm = TRUE)
  )
comparison_data <- rbind(
  comparison_data,
  temp
)
temp <- 
  data.frame(
    comparison = "IQR",
    RF = apply(boot_resids1,1, IQR,na.rm = TRUE),
    LR = apply(lm_bootresids_OOB,1, IQR,na.rm = TRUE)
  )
comparison_data <- rbind(
  comparison_data,
  temp
)
ggplot(data = comparison_data) + geom_point(aes(x  = abs(RF), y = abs(LR))) + facet_wrap(~comparison,scales = "free") +
  geom_abline(slope = 1,intercept = 0)
}

{
  bootstrap_residuals <- function(orig_data, caret_fit, n_resamples = 1000, strata = NULL, 
                                  ncores = 1, verbose = FALSE, force_OOB = TRUE) 
  {
    train.formula <- caret:::train.formula
    resample_func <- function(x) {
      if (is.null(strata)) {
        bs_samples <- sample(1:nrow(orig_data), replace = TRUE)
      }
      else {
        groups <- sort(unique(strata))
        bs_samples <- list()
        for (g in groups) {
          group <- which(strata == g)
          bs_samples[[g]] <- sample(group, replace = TRUE)
        }
        bs_samples <- unlist(bs_samples)
      }
      formula <- caret_fit$terms
      response <- caret_fit$terms[[2]]
      bs_data <- orig_data[bs_samples, ]
      caret_call <- caret_fit$call
      new_call <- pryr::modify_call(caret_call, new_args = list(data = quote(bs_data), 
                                                                trControl = quote(trainControl(method = "none")), 
                                                                tuneGrid = quote(caret_fit$bestTune)))
      bs_caret_fit <- eval(new_call)
      pred <- predict(bs_caret_fit, newdata = orig_data)
      if (force_OOB) {
        in_sample <- which(1:nrow(orig_data) %in% bs_samples)
        pred[in_sample] <- NA
      }
      resid <- orig_data[[response]] - pred
      resid
    }
    t1 <- proc.time()
    myCluster <- parallel::makeCluster(ncores)
    parallel::clusterCall(cl = myCluster, fun = function(x) {
      lapply(X = c("caret", "pryr"), FUN = function(x) library(x, 
                                                               character.only = TRUE))
    })
    parallel::clusterExport(cl = myCluster, varlist = c("orig_data", 
                                                        "caret_fit", "train.formula", "strata"), envir = environment())
    resid <- cluster_lapply2(x = 1:n_resamples, func = resample_func, 
                             cluster = myCluster)
    parallel::stopCluster(myCluster)
    t2 <- proc.time()
    if (verbose) {
      message("Completed ", n_resamples, " bootstrap resamples in ", 
              round((t2 - t1)[3]/60, 2), " minutes.")
    }
    resid_save.mat <- do.call(what = "cbind", resid)
    resid_save.mat
  }
  
  
}

{
  lm_bootresids_OOB <- bootstrap_residuals(
    orig_data = data_list[[r]],
    caret_fit = caret_lm_fit,
    n_resamples = 3000,
    ncores = 7,
    verbose = TRUE,
    force_OOB = TRUE,
    strata = clusters
  )
  
  boot_resids1 <- bootstrap_residuals(orig_data = data_list[[r]],
                                      caret_fit = caret_fit,
                                      n_resamples = 3000,
                                      ncores = 7,
                                      verbose = TRUE,
                                      force_OOB = TRUE,
                                      strata = clusters
  )
  
  rankings <- get_bootstrap_rankings(
    boot_resids = boot_resids1,
    alpha_level = 0.9,
    group = data_list[[r]]$cluster,n_samples = 50000)
  
  lm_rankings <- get_bootstrap_rankings(
    boot_resids = lm_bootresids_OOB,
    alpha_level = 0.9,
    group = data_list[[r]]$cluster,n_samples = 50000)
  
  rankings$residual_stats$method = "RF"
  lm_rankings$residual_stats$method = "LR"
  
  full_ranks <- rbind(rankings$residual_stats,lm_rankings$residual_stats)
  
  
  full_ranks %>% select(method,mean_rank) %>% 
    mutate(id = rep(1:200,times = 2)) %>% 
    pivot_wider(id_cols = "id",names_from = method,values_from = mean_rank  ) %>%
    ggplot(aes(x=RF, y=LR)) + geom_point() + geom_abline(slope=1,intercept =0)
  
  
  full_ranks %>% select(method,lower_rank,upper_rank) %>% 
    mutate(id = rep(1:200,times = 2),
           value = upper_rank - lower_rank) %>% 
    pivot_wider(id_cols = "id",names_from = method,values_from = value) %>%
    ggplot(aes(x=RF, y=LR)) + geom_point() + geom_abline(slope=1,intercept =0)
  
  mean(abs(apply(boot_resids1,1,mean,na.rm = TRUE)))
  mean(abs(apply(lm_bootresids_OOB,1,mean,na.rm = TRUE)))
  
  blah <- 
    full_ranks %>% select(method,lower_rank,upper_rank) %>% 
    mutate(id = rep(1:200,times = 2),
           value = upper_rank - lower_rank) %>% 
    pivot_wider(id_cols = "id",names_from = method,values_from = value)
  mean(blah$RF - blah$LR < 0)
  hist(blah$RF - blah$LR)
  median(blah$RF - blah$LR)
  
  which(blah$RF)
  
  which.max(blah$LR)
  
}

{
  # plot(predict(caret_fit$finalModel)-data_list[[r]][[r]])
  par(mfrow = c(4,4))
  for(i in 2:ncol(data_list[[r]])){
    plot(y= predict(caret_fit$finalModel)-data_list[[r]][[r]],x = data_list[[r]][[i]],main = colnames(data_list[[r]])[i])
    abline(0,0)
    lines(predict)
  }
  
  
  
  i = 1
  # plot(predict(caret_fit$finalModel)-data_list[[r]][[r]])
  par(mfrow = c(4,4))
  for(i in 2:ncol(data_list[[r]])){
    plot(y= predict(caret_lm_fit)-data_list[[r]][[r]],x = data_list[[r]][[i]],main = colnames(data_list[[r]])[i])
    abline(0,0)
    
  }
  # plot(y= predict(caret_lm_fit)-data_list[[r]][[r]],x = data_list[[r]][[i]])
  
  predict_data <- as.data.frame(Benchtools:::newdata_caret_prepare(data_list[[r]],train_fit = caret_lm_fit))
  predict_data$LR_residual = predict(caret_lm_fit)-data_list[[r]][[r]]
  predict_data$RF_residual = predict(caret_fit$finalModel)-data_list[[r]][[r]]
  predict_data <- predict_data %>% reshape2::melt(id= c("LR_residual","RF_residual"))
  important_variables <- c("covariate8","covariate26","covariate14","covariate10")
  p1 <- ggplot(data = predict_data %>% filter(variable %in% important_variables),mapping = aes(x = value, y = LR_residual)) + geom_point() + 
    geom_hline(yintercept = 0) +
    geom_smooth(method = "gam") + facet_wrap(~variable,scales = "free",nrow = 1) + theme_bw() + labs(y = "Residual (Linear Regression)",x = "Covariate value")
  
  p2 <- ggplot(data = predict_data%>% filter(variable %in% important_variables),mapping = aes(x = value, y = RF_residual)) + geom_point() + 
    geom_hline(yintercept = 0) +
    geom_smooth(method = "gam") + facet_wrap(~variable,scales = "free",nrow = 1) + theme_bw() + labs(y = "Residual (Random Forest)",x = "Covariate value")
  
  pdf(file = "slides/figures/residual-plots.pdf",width = 6,height = 3,onefile = FALSE)
  ggpubr::ggarrange(p1,p2,nrow = 2,ncol = 1)
  dev.off()
  
  pdf(file = "benchmarking/figs/residual-plots.pdf",width = 8,height = 5,onefile = FALSE)
  ggpubr::ggarrange(p1,p2,nrow = 2,ncol = 1)
  dev.off()
  predict_data %>% pivot_longer(LR_residual,RF_residual,names_to = "model",values_to = "value")
}