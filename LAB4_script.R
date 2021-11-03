
# 1.Read file

# 2.Check validity of data

# 3.Add data to summary table

# 4.Preliminary visualization

# 5.Define ensemble of models


# 6.1 Obtain initial parameters with linear regression


# 6.2 Perform nonlinear regression

# 6.0 Check homoscesdasticity holds

# 6.3 Get AIC and other necessary coefs

# 6.4 Build results tables

# 7. Plot real data vs best fit 

# Final script

langs_summary <- function(languages){
  table_1 <- data.frame(language=languages,
                        N = numeric(length(languages)),
                        mu_n = numeric(length(languages)),
                        sd_n = numeric(length(languages)),
                        mu_k = numeric(length(languages)),
                        sd_k = numeric(length(languages)))
  
  for (l in 1:length(languages)) {
    lang <- languages[l]
    path <- paste("./data/",paste(lang,"_dependency_tree_metrics.txt",sep = ""),sep = "")
    raw_data <- read.table(path,header = FALSE,sep = " ",fill = TRUE,quote="", comment.char="") 
    colnames(raw_data) = c("vertices","degree_2nd_moment", "mean_length")
    raw_data = raw_data[order(raw_data$vertices),]
    # check validity
    valid_data <- raw_data %>% mutate(check = case_when(        
      (degree_2nd_moment >= (4-6)/vertices & degree_2nd_moment <= vertices-1) ~ "valid",
      TRUE ~ "invalid")) %>% filter(check == "valid")
    
    invalids = nrow(raw_data)-nrow(valid_data)
    print(paste("Invalid values for",lang,":",invalids))
    
    table_1$N[l] <- nrow(valid_data)
    table_1$mu_n[l] <- mean(valid_data$vertices)
    table_1$sd_n[l] <- sd(valid_data$vertices)
    table_1$mu_k[l] <- mean(valid_data$degree_2nd_moment)
    table_1$sd_k[l] <- sd(valid_data$degree_2nd_moment)
  }
  
  table_1
}

ensemble_of_models <- function(language){
  # Read file
  file <- paste("./data/", language, 
                "_dependency_tree_metrics.txt",
                sep = "")
  lang_data = read.table(file, header = FALSE)
  colnames(lang_data) = c("vertices","degree_2nd_moment", "mean_length")
  lang_data = lang_data[order(lang_data$vertices), ]
  
  # Obtain mean of lang
  mean_lang = aggregate(lang_data, list(lang_data$vertices), mean)
  
  # MODEL 0
  ## Get AIC, standard error
  n <- length(mean_lang$vertices)
  RSS <- sum((mean_lang$degree_2nd_moment - (1-1/n)*(5-6/n)))
  p <- 0
  s_0 <- sqrt(RSS/(n - p))
  AIC_0 <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)
  
  # MODEL 1
  ## Obtain initial parameters with linear regression
  linear_model_1 = lm(log(degree_2nd_moment)~log(vertices), mean_lang)
  b_initial = coef(linear_model_1)[2]
  
  ## Perform nonlinear regression
  nonlinear_model_1 = nls(degree_2nd_moment~(vertices/2)^b,data=mean_lang,
                        start = list(b = b_initial), trace = FALSE)
  
  ## Get AIC, RSS and obtained coefficients
  s_1 <- sqrt(deviance(nonlinear_model_1)/df.residual(nonlinear_model_1))
  AIC_1 <- AIC(nonlinear_model_1)
  b_1 <- coef(nonlinear_model_1)["b"]
  
  # MODEL 2
  ## Obtain initial parameters with linear regression
  linear_model_2 = lm(log(degree_2nd_moment)~log(vertices), mean_lang)
  a_initial = exp(coef(linear_model_2)[1])
  b_initial = coef(linear_model_2)[2]
  
  
  ## Perform nonlinear regression
  nonlinear_model_2 = nls(degree_2nd_moment~a*vertices^b,data=mean_lang,
                          start = list(a = a_initial, b = b_initial), trace = FALSE)
  
  ## Get AIC, RSS and obtained coefficients
  s_2 <- sqrt(deviance(nonlinear_model_2)/df.residual(nonlinear_model_2))
  AIC_2 <- AIC(nonlinear_model_2)
  a_2 <- coef(nonlinear_model_2)["a"]
  b_2 <- coef(nonlinear_model_2)["b"]
  
  # MODEL 3
  ## Obtain initial parameters with linear regression
  linear_model_3 = lm(log(degree_2nd_moment)~vertices, mean_lang)
  a_initial = exp(coef(linear_model_3)[1])
  c_initial = coef(linear_model_3)[2]
  
  
  ## Perform nonlinear regression
  nonlinear_model_3 = nls(degree_2nd_moment~a*exp(c*vertices),data=mean_lang,
                          start = list(a = a_initial, c= c_initial), trace = FALSE)
  
  ## Get AIC, RSS and obtained coefficients
  s_3 <- sqrt(deviance(nonlinear_model_3)/df.residual(nonlinear_model_3))
  AIC_3 <- AIC(nonlinear_model_3)
  a_3 <- coef(nonlinear_model_3)["a"]
  c_3 <- coef(nonlinear_model_3)["c"]
  
  
  # Build final results vectors
  coefs <- c(language, 
              b_1,
              a_2, b_2,
              a_3, c_3)
  aics <- c(language, AIC_0, AIC_1, AIC_2, AIC_3)
  
  best_aic <- min(c(AIC_0, AIC_1, AIC_2, AIC_3))
  
  aics_diff <- c(language, abs(AIC_0-best_aic),
                 abs(AIC_1-best_aic), abs(AIC_2-best_aic), 
                 abs(AIC_3-best_aic))
  
  std_error <- c(language, s_0, s_1, s_2, s_3)
  
  return(list("coefficients" = coefs ,"aics"=aics ,"aics_diff" = aics_diff, 
              "std_error" = std_error))
}


languages <- c("Arabic", "Basque", "Catalan", "Chinese", "Czech", "English", "Greek", "Hungarian", "Italian", "Turkish")

summary_table <- langs_summary(languages)

data <- ensemble_of_models("Catalan")
