
# read a language file
read_language <- function(lang) {
  
  path <- paste("./data/",paste(lang,"_dependency_tree_metrics.txt",sep = ""),sep = "")
  raw_data <- read.table(path,header = FALSE,sep = " ",fill = TRUE,quote="", comment.char="") 
  colnames(raw_data) = c("vertices","degree_2nd_moment", "mean_length")
  raw_data = raw_data[order(raw_data$vertices),]
  
  # check validity
  valid_data <- raw_data %>% mutate(check = case_when(        
    (degree_2nd_moment >= (4-6)/vertices & degree_2nd_moment <= vertices-1) ~ "valid",
    TRUE ~ "invalid")) %>% filter(check == "valid") %>% select(-check)
  return(valid_data)
}


# plot data with line at mean+4*sd
# return outlier free dataframe if remove=TRUE
find_outliers <- function(data,lang,remove,plot) {
  
  mean <- mean(data$degree_2nd_moment)
  sd <- sd(data$degree_2nd_moment)
  
  if (plot == TRUE) {
    plot(data$vertices, data$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment")
    abline(h=mean+4*sd, col='green')  
  }

  
  # count
  n_outliers <- length(which(data$degree_2nd_moment > mean+4*sd))
  print(paste("outliers removed for",lang,":",n_outliers))
  
  if (remove == TRUE) {
    # remove
    data <- data %>% mutate(outlier = case_when(
      degree_2nd_moment > mean+4*sd ~ "out",
      TRUE ~ "in" )) %>% filter(outlier == "in") %>% select(-outlier)
    
    return(data)
  }
  
}




library(minpack.lm)
library(sandwich)

nlsvcovhc <- function(nlslmout) {
  
  # uses output of nlsLM() of the minpack.lm package to get an asymptotic 
  # covariance matrix without assuming homoscedasticity
  
  b <- coef(nlslmout)
  m <- nlslmout$m
  resid <- m$resid()
  hmat <- m$gradient()
  fakex <- hmat
  fakey <- resid + hmat %*% b
  lmout <- lm(fakey ~ fakex - 1)
  vcovHC(lmout)
}


get_coeff_CI <- function(robust) {
  c_CI_list <- c()
  for (i in 1:length(languages)) {
    lang <- languages[i]
    data_raw <- read_language(lang)
    data <- find_outliers(data_raw,lang, remove=T,plot=F)
    data <- aggregate(data, list(data$vertices), mean)
    n <- nrow(data)
    
    linear_model_3 = lm(log(degree_2nd_moment)~vertices, data)
    a_initial = exp(coef(linear_model_3)[1])
    c_initial = coef(linear_model_3)[2]
    d_initial <- 0
    
    if (robust == TRUE) {
      nonlinear_model_3p = nlsLM(degree_2nd_moment~a*exp(c*vertices)+d,data=data,
                               start = list(a = a_initial, c= c_initial,
                                            d = d_initial), 
                               trace = FALSE,
                               control = list(maxiter=1000),
                               algorithm = 'port',
                               lower=c(0,0,-20)
                               )
      sd <- diag(sqrt(nlsvcovhc(nonlinear_model_3p)))[2]
      c_CI_list[i] <-  1.96*(sd/sqrt(n))
      
      
    }
    
    else {
      nonlinear_model_3p = nls(degree_2nd_moment~a*exp(c*vertices)+d,data=data,
                               start = list(a = a_initial, c= c_initial,
                                            d = d_initial), 
                               trace = FALSE,
                               control = list(maxiter=1000),
                               algorithm = 'port',
                               lower=c(0,0,-20)
                                )
      sd <- diag(sqrt(vcov(nonlinear_model_3p)))[2]
      c_CI_list[i] <-  1.96*(sd/sqrt(n))
    }
    
  }
  return(data.frame("half_CI_c"=c_CI_list))
}





