
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
find_outliers <- function(data,remove) {
  
  mean <- mean(data$degree_2nd_moment)
  sd <- sd(data$degree_2nd_moment)
  
  plot(data$vertices, data$degree_2nd_moment, xlab = "vertices", ylab = "degree 2nd moment")
  abline(h=mean+4*sd, col='green')
  
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


# plot model with confidence intervals
plot_with_CI <- function(data, nls_model, conf_level) {
  plotFit(nls_model, interval = "confidence", level = conf_level, data, adjust= 'none',
          shade = TRUE, col.conf = "pink", col.fit = "red")
  lines(data$vertices,4-6/data$vertices, col = "blue")
  lines(data$vertices,data$vertices-1, col = "blue")
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


