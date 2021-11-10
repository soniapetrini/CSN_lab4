
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

source("funs.R")

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

preliminary_plots <- function(){
  # Read file
  file <- paste("./data/", language, 
                "_dependency_tree_metrics.txt",
                sep = "")
  lang_data = read.table(file, header = FALSE)
  colnames(lang_data) = c("vertices","degree_2nd_moment", "mean_length")
  lang_data = lang_data[order(lang_data$vertices), ]
  
  
  plot(lang_data$vertices, lang_data$degree_2nd_moment,
       xlab = "vertices", ylab = "degree_2nd_moment")
  
  plot(log(lang_data$vertices), log(lang_data$mean_length),
         xlab = "log(vertices)", ylab = "log(degree_2nd_moment)")
  
  mean_lang = aggregate(lang_data, list(lang_data$vertices), mean)
  
  plot(mean_lang$vertices, mean_lang$degree_2nd_moment,
       xlab = "vertices", ylab = "mean degree_2nd_moment")

  plot(log(mean_lang$vertices), log(mean_lang$degree_2nd_moment),
         xlab = "log(vertices)", ylab = "log(mean degree_2nd_moment)")
  
  plot(log(lang_data$vertices), log(lang_data$degree_2nd_moment),
       xlab = "vertices", ylab = "degree_2nd_moment")
  lines(log(mean_lang$vertices),log(mean_lang$degree_2nd_moment), col = "green")
  # lines(log(mean_lang$vertices),log((mean_lang$vertices+1)/3), col = "red")
  # cambiar por formula de k^2
  
  
  plot(lang_data$vertices, lang_data$degree_2nd_moment,
       xlab = "vertices", ylab = "degree 2nd moment")
  lines(mean_lang$vertices,mean_lang$degree_2nd_moment, col = "green")
  lines(lang_data$vertices,
        (1 - 1/lang_data$vertices)*(5 - 6/lang_data$vertices), col = "red")
  lines(lang_data$vertices,4-6/lang_data$vertices, col = "blue")
  lines(lang_data$vertices,lang_data$vertices-1, col = "blue")
  
}

ensemble_of_models <- function(language){
  # Read file
  file <- paste("./data/", language, 
                "_dependency_tree_metrics.txt",
                sep = "")
  lang_data = read.table(file, header = FALSE)
  colnames(lang_data) = c("vertices","degree_2nd_moment", "mean_length")
  lang_data = lang_data[order(lang_data$vertices), ]
  
  # find and remove outliers
  #lang_data <- find_outliers(lang_data,remove=TRUE)
  
  # Obtain mean of lang
  mean_lang = aggregate(lang_data, list(lang_data$vertices), mean)
  
  # MODEL 0
  ## Get AIC, standard error
  RSS <- sum((mean_lang$degree_2nd_moment - (1-1/mean_lang$vertices)*(5-6/mean_lang$vertices))^2)
  n <- length(mean_lang$vertices)
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
  
  # MODEL 1+
  ## Perform nonlinear regression
  # Here we use the previous initial b obtained with linear regression and a
  # random starting value for d
  d_initial <- 7
  model_1plus = nls(degree_2nd_moment~(vertices/2)^b + d, data=mean_lang, 
                    start = list(b = b_initial, d = d_initial))
  ## Get AIC, RSS and obtained coefficients
  b_1p <- coef(model_1plus)["b"]
  d_1p <- coef(model_1plus)["d"]
  s_1p <- sqrt(deviance(model_1plus)/df.residual(model_1plus))
  AIC_1p <- AIC(model_1plus)
  
  
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
  
  # MODEL 2+
  ## Perform nonlinear regression
  # Here we use the previous initial a,b obtained with linear regression and a
  # random starting value for d
  d_initial <- 0
  model_2plus = nls(degree_2nd_moment~a*(vertices)^b + d, data=mean_lang, 
                    start = list(a = a_initial, b = b_initial, d = d_initial),
                    control = list(maxiter=1000),
                    algorithm = 'port',
                    lower=c(0,0,-20),
                    upper=c(20,5,20))
  ## Get AIC, RSS and obtained coefficients
  a_2p <- coef(model_2plus)["a"]
  b_2p <- coef(model_2plus)["b"]
  d_2p <- coef(model_2plus)["d"]
  s_2p <- sqrt(deviance(model_2plus)/df.residual(model_2plus))
  AIC_2p <- AIC(model_2plus)
  
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
  
  
  # MODEL 3+
  ## Perform nonlinear regression
  d_initial <- 0
  nonlinear_model_3p = nls(degree_2nd_moment~a*exp(c*vertices)+d,data=mean_lang,
                          start = list(a = a_initial, c= c_initial,
                                       d = d_initial), 
                          trace = FALSE,
                          control = list(maxiter=1000),
                          algorithm = 'port',
                          lower=c(0,0,-20)
                          )
  
  ## Get AIC, RSS and obtained coefficients
  s_3p <- sqrt(deviance(nonlinear_model_3p)/df.residual(nonlinear_model_3p))
  AIC_3p <- AIC(nonlinear_model_3p)
  a_3p <- coef(nonlinear_model_3p)["a"]
  c_3p <- coef(nonlinear_model_3p)["c"]
  d_3p <- coef(nonlinear_model_3p)["d"]
  
  
  # MODEL 4
  ## Obtain initial parameters with linear regression
  linear_model_4 = lm(log(degree_2nd_moment)~log(vertices), mean_lang)
  a_initial = coef(linear_model_3)[1]
  
  
  ## Perform nonlinear regression
  nonlinear_model_4 = nls(degree_2nd_moment~a*log(vertices),
                          data=mean_lang,
                          start = list(a = a_initial), 
                          trace = FALSE)
  
  ## Get AIC, RSS and obtained coefficients
  s_4 <- sqrt(deviance(nonlinear_model_4)/df.residual(nonlinear_model_4))
  AIC_4 <- AIC(nonlinear_model_4)
  a_4 <- coef(nonlinear_model_4)["a"]
  
  
  # MODEL 4+
  d_initial <- 0
  ## Perform nonlinear regression
  nonlinear_model_4p = nls(degree_2nd_moment~a*log(vertices) + d,
                          data=mean_lang,
                          start = list(a = a_initial, d=d_initial), 
                          trace = FALSE)
  
  ## Get AIC, RSS and obtained coefficients
  s_4p <- sqrt(deviance(nonlinear_model_4p)/df.residual(nonlinear_model_4p))
  AIC_4p <- AIC(nonlinear_model_4p)
  a_4p <- coef(nonlinear_model_4p)["a"]
  d_4p <- coef(nonlinear_model_4p)["d"]
  
  # Build final results vectors
  coefs <- c(language, 
              round(b_1,3),
              round(a_2,3), round(b_2,3),
              round(a_3,3), round(c_3,3),
              round(a_4,3),
              round(b_1p,3), round(d_1p,3),
              round(a_2p,3), round(b_2p,3), round(d_2p,3),
              round(a_3p,3), round(c_3p,3), round(d_3p,3),
              round(a_4p,3), round(d_4p,3))
  aics <- c(language, round(AIC_0,3), round(AIC_1,3), round(AIC_2,3), round(AIC_3,3), round(AIC_4,3), 
            round(AIC_1p,3), round(AIC_2p,3), round(AIC_3p,3), round(AIC_4p,3))
  
  best_aic <- min(c(AIC_0, AIC_1, AIC_2, AIC_3, AIC_4, 
                    AIC_1p, AIC_2p, AIC_3p, AIC_4p))
  
  aics_diff <- c(language, round(abs(AIC_0-best_aic),3),
                 round(abs(AIC_1-best_aic),3),
                 round(abs(AIC_2-best_aic),3), 
                 round(abs(AIC_3-best_aic),3),
                 round(abs(AIC_4-best_aic),3),
                 round(abs(AIC_1p-best_aic),3),
                 round(abs(AIC_2p-best_aic),3),
                 round(abs(AIC_3p-best_aic),3),
                 round(abs(AIC_4p-best_aic),3))
  
  std_error <- c(language, round(s_0,3), round(s_1,3), round(s_2,3), round(s_3,3), round(s_4,3),
                 round(s_1p,3), round(s_2p,3), round(s_3p,3), round(s_4p,3))
  
  return(list("coefficients" = coefs ,"aics"=aics ,"aics_diff" = aics_diff, 
              "std_error" = std_error))
}

# Initialize constants
languages <- c("Arabic", "Basque", "Catalan", "Chinese", "Czech", "English", "Greek", "Hungarian", "Italian", "Turkish")

# Initialize dataframes
coefficients_df <- data.frame(matrix(ncol = 17, nrow = 0))
coefficients_df_names <- c('Language', '1: b', '2: a', '2: b', '3: a',
                     '3: c', '4: a', '1+: b', '1+: d', '2+: a',
                     '2+: b', '2+: d', '3+: a', '3+: c', '3+: d',
                     '4+: a', '4+: d')
colnames(coefficients_df) <- coefficients_df_names

residual_err_df <- data.frame(matrix(ncol = 10, nrow = 0))
residual_err_df_names <- c('Language', '0', '1', '2', '3', '4', '1+',
                           '2+', '3+', '4+')
colnames(residual_err_df) <- residual_err_df_names

aic_df <- data.frame(matrix(ncol = 10, nrow = 0))
aic_df_names <- c('Language', '0', '1', '2', '3', '4', '1+',
                           '2+', '3+', '4+')
colnames(aic_df) <- aic_df_names

aic_diff_df <- data.frame(matrix(ncol = 10, nrow = 0))
aic_diff_df_names <- c('Language', '0', '1', '2', '3', '4', '1+',
                  '2+', '3+', '4+')
colnames(aic_diff_df) <- aic_diff_df_names

# Show summary
summary_table <- langs_summary(languages)

# Show model fitting results
for (language in languages){
  lang_data <- ensemble_of_models(language)

  # Add models params to df
  language_coeffs_df = data.frame(t(unlist(lang_data$coefficients)))
  names(language_coeffs_df) = coefficients_df_names
  coefficients_df <- rbind(coefficients_df,language_coeffs_df)
  
  # Add models residual errors to df
  lang_err_df = data.frame(t(unlist(lang_data$std_error)))
  names(lang_err_df) = residual_err_df_names
  residual_err_df <- rbind(residual_err_df,lang_err_df)
  
  # Add models AIC to df
  lang_aic_df = data.frame(t(unlist(lang_data$aics)))
  names(lang_aic_df) = aic_df_names
  aic_df <- rbind(aic_df,lang_aic_df)
  
  # Add models AIC diff to df
  lang_aic_diff_df = data.frame(t(unlist(lang_data$aics_diff)))
  names(lang_aic_diff_df) = aic_diff_df_names
  aic_diff_df <- rbind(aic_diff_df,lang_aic_diff_df)
}

coefficients_df

residual_err_df

aic_df

aic_diff_df

