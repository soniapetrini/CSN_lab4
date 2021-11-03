library(dplyr)

degree_seq_summary_df <- data.frame(matrix(ncol = 6, nrow = 0))
degree_seq_summary_df_names <-  c('Language', 'N', 'n_mean', 'n_std', 'x_mean', 
                                  'x_std')
colnames(degree_seq_summary_df) <- degree_seq_summary_df_names

language = "Catalan"

# 1.Read file
file <- paste("./data/", language, 
              "_dependency_tree_metrics.txt",
              sep = "")
lang_data = read.table(file, header = FALSE)
colnames(lang_data) = c("vertices","degree_2nd_moment", "mean_length")
lang_data = lang_data[order(lang_data$vertices), ]
# 2.Check validity of data

# 3.Add data to summary table
N <- nrow(lang_data)
n_mean <- mean(lang_data$vertices)
n_std <- sd(lang_data$vertices)
x_mean <- mean(lang_data$degree_2nd_moment)
x_std <- sd(lang_data$degree_2nd_moment)

language_summary_df = data.frame(t(unlist(c(language,N,n_mean,n_std,x_mean,x_std))))
names(language_summary_df) = degree_seq_summary_df_names
degree_seq_summary_df <- rbind(degree_seq_summary_df,language_summary_df)

# 4.Preliminary visualization

# 5.Define ensemble of models

##### Model 1
# 6.1 Obtain initial parameters with linear regression
linear_model = lm(log(degree_2nd_moment)~log(vertices/2), lang_data)
a_initial = exp(coef(linear_model)[1])
b_initial = coef(linear_model)[2]

# 6.2 Perform nonlinear regression
nonlinear_model = nls(degree_2nd_moment~a*vertices^b,data=lang_data,
                      start = list(a = a_initial, b = b_initial), trace = TRUE)
# 6.0 Check homoscesdasticity holds

# 6.3 Get AIC and other necessary coefs
RSS_2 <- deviance(nonlinear_model)
AIC_2 <- AIC(nonlinear_model)

a_2 <- coef(nonlinear_model)["a"]
b_2 <- coef(nonlinear_model)["b"]



##### Model 2
# 6.1 Obtain initial parameters with linear regression
linear_model = lm(log(degree_2nd_moment)~log(vertices), lang_data)
a_initial = exp(coef(linear_model)[1])
b_initial = coef(linear_model)[2]

# 6.2 Perform nonlinear regression
nonlinear_model = nls(degree_2nd_moment~a*vertices^b,data=lang_data,
                      start = list(a = a_initial, b = b_initial), trace = TRUE)
# 6.0 Check homoscesdasticity holds

# 6.3 Get AIC and other necessary coefs
RSS_2 <- deviance(nonlinear_model)
AIC_2 <- AIC(nonlinear_model)

a_2 <- coef(nonlinear_model)["a"]
b_2 <- coef(nonlinear_model)["b"]

# 6.4 Build results tables

# 7. Plot real data vs best fit 