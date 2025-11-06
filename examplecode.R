# Load required packages
library(survival)
library(flexsurv)
library(dplyr)
library(ggplot2)


# Load the data
methadone <- read.csv("methadone.csv")

# (a) Compute average follow-up time and proportion censored
# Note: 'event' is assumed to be 1 = event occurred, 0 = censored

# Average follow-up time
avg_followup <- mean(methadone$time)
cat("Average follow-up time (days):", avg_followup, "\n")

# Proportion of censored observations
prop_censored <- mean(methadone$event == 0)
cat("Proportion censored:", prop_censored, "\n")

# (b): Fit Exponential, Weibull, and Generalized Gamma models and Report parameter estimates and log-likelihoods
# Create survival object
surv_obj <- Surv(time = methadone$time, event = methadone$event)

# Fit models
fit_exp <- flexsurvreg(surv_obj ~ 1, dist = "exponential")
fit_weibull <- flexsurvreg(surv_obj ~ 1, dist = "exponential")
fit_gengamma <- flexsurvreg(surv_obj ~ 1, dist = "gengamma")


# (c) Plots
# Create a grid of time values
time_grid <- seq(0, max(methadone$time), length.out = 200)

# Get predicted survival from each model
pred_exp <- summary(fit_exp, t = time_grid, type = "survival")[[1]]
pred_weibull <- summary(fit_weibull, t = time_grid, type = "survival")[[1]]
pred_gengamma <- summary(fit_gengamma, t = time_grid, type = "survival")[[1]]

# Combine into one data frame
df_models <- bind_rows(
  data.frame(time = pred_exp$time, surv = pred_exp$est, model = "Exponential"),
  data.frame(time = pred_gengamma$time, surv = pred_gengamma$est, model = "Generalized Gamma")
)

# Get Kaplan-Meier survival estimate
km_fit <- survfit(surv_obj ~ 1)
df_km <- data.frame(time = km_fit$time, surv = km_fit$surv, model = "Kaplan-Meier")

# Plot all together
ggplot() +
  geom_step(data = df_km, aes(x = time, y = surv, color = model), size = 1) +
  geom_line(data = df_models, aes(x = time, y = surv, color = model), size = 1) +
  labs(title = "Survival Function: Parametric vs Nonparametric Estimates",
       x = "Time (days)",
       y = "Survival Probability") +
  theme_minimal() +
  scale_color_manual(values = c(
    "Kaplan-Meier" = "black",
    "Exponential" = "red",
    "Weibull" = "blue",
    "Generalized Gamma" = "darkgreen"
  )) +
  theme(legend.title = element_blank())


# (d) Likelihood Ratio Test: Is Weibull a simplification of GG?
# Extract log-likelihoods
loglik_weibull <- fit_weibull$loglik
loglik_gengamma <- fit_gengamma$loglik

# Compute likelihood ratio statistic
lr_stat <- 2 * (loglik_gengamma - loglik_weibull)

# Degrees of freedom: difference in number of parameters
df <- fit_gengamma$npars - fit_weibull$npars

# Compute p-value
p_value <- pchisq(lr_stat, df = df, lower.tail = FALSE)

# Print results
cat("Likelihood Ratio Statistic:", round(lr_stat, 3), "\n")
cat("Degrees of Freedom:", df, "\n")
cat("P-value:", round(p_value, 4), "\n")

# Interpretation
if(p_value < 0.05){
  cat("Conclusion: The generalized gamma model fits significantly better. The Weibull model is NOT an adequate simplification.\n")
} else {
  cat("Conclusion: No significant improvement with generalized gamma. The Weibull model is an adequate simplification.\n")
}


# (e) estimate and 95%CI of median time to exit and probability of exit by 1 and 2 years
### (i) Median time to exit from maintenance
cat("\n(i) Median survival time:\n")
fitparametric(surv_obj, dist = "weibull", feature = "quantile", pi = 0.5)

### (ii) Probability of survival beyond 1 year (365 days)
cat("\n(ii) Probability of survival beyond 1 year:\n")
fitparametric(surv_obj, dist = "weibull", feature = "survival", t = 365)

### (iii) Conditional survival: survive to 2 years given survived 1 year
cat("\n(iii) Conditional survival P(T > 730 | T > 365):\n")
fitparametric(surv_obj, dist = "weibull", feature = "condsurvival", t = 730, t0 = 365)



# (f) Likelihood Ratio Test: Is Exponential a simplification of Weibull?
# Extract log-likelihoods
loglik_exp <- fit_exp$loglik
loglik_weibull <- fit_weibull$loglik

# Compute likelihood ratio statistic
lr_stat_exp_weib <- 2 * (loglik_weibull - loglik_exp)

# Degrees of freedom (difference in number of parameters)
df <- fit_weibull$npars - fit_exp$npars

# Compute p-value
p_value_exp_weib <- pchisq(lr_stat_exp_weib, df = df, lower.tail = FALSE)

# Print results
cat("Likelihood Ratio Statistic:", round(lr_stat_exp_weib, 3), "\n")
cat("Degrees of Freedom:", df, "\n")
cat("P-value:", round(p_value_exp_weib, 4), "\n")

# Interpretation
if(p_value_exp_weib < 0.05){
  cat("Conclusion: The Weibull model fits significantly better. The Exponential model is NOT an adequate simplification.\n")
} else {
  cat("Conclusion: No significant improvement with Weibull. The Exponential model is an adequate simplification.\n")
}



# (g) Compare Clinic 1 vs Clinic 2 using exponential models
# Subset data
clinic1_data <- methadone %>% filter(clinic == 1)
clinic2_data <- methadone %>% filter(clinic == 2)

# Create survival objects
surv_clinic1 <- Surv(time = clinic1_data$time, event = clinic1_data$event)
surv_clinic2 <- Surv(time = clinic2_data$time, event = clinic2_data$event)

# Fit exponential models using fitparametric
fit1 <- fitparametric(surv_clinic1, dist = "exp")
fit2 <- fitparametric(surv_clinic2, dist = "exp")

# Model without group (null)
surv_all <- Surv(time = methadone$time, event = methadone$event)
fit_null <- flexsurvreg(Surv(time, event) ~ 1, data = methadone, dist = "exponential")

# Model with group (clinic as factor)
fit_alt <- flexsurvreg(Surv(time, event) ~ factor(clinic), data = methadone, dist = "exponential")

loglik_null <- fit_null$loglik
loglik_alt <- fit_alt$loglik

lrt_stat <- 2 * (loglik_alt - loglik_null)
df <- fit_alt$npars - fit_null$npars
p_val <- pchisq(lrt_stat, df = df, lower.tail = FALSE)

cat("\nLikelihood Ratio Test for clinic effect:\n")
cat("LRT statistic:", round(lrt_stat, 3), "\n")
cat("Degrees of freedom:", df, "\n")
cat("P-value:", round(p_val, 4), "\n")


# (h)
# Subset data
prison0_data <- methadone %>% filter(prison == 0)
prison1_data <- methadone %>% filter(prison == 1)

# Survival objects
surv_prison0 <- Surv(time = prison0_data$time, event = prison0_data$event)
surv_prison1 <- Surv(time = prison1_data$time, event = prison1_data$event)

# Fit exponential models using fitparametric
fit_p0 <- fitparametric(surv_prison0, dist = "exp")
fit_p1 <- fitparametric(surv_prison1, dist = "exp")

# Fit models with and without prison history
fit_null_prison <- flexsurvreg(Surv(time, event) ~ 1, data = methadone, dist = "exponential")
fit_alt_prison <- flexsurvreg(Surv(time, event) ~ prison, data = methadone, dist = "exponential")

# LRT
loglik_null_p <- fit_null_prison$loglik
loglik_alt_p <- fit_alt_prison$loglik

lrt_stat_p <- 2 * (loglik_alt_p - loglik_null_p)
df_p <- fit_alt_prison$npars - fit_null_prison$npars
p_val_p <- pchisq(lrt_stat_p, df = df_p, lower.tail = FALSE)

cat("\nLikelihood Ratio Test for prison history effect:\n")
cat("LRT statistic:", round(lrt_stat_p, 3), "\n")
cat("Degrees of freedom:", df_p, "\n")
cat("P-value:", round(p_val_p, 4), "\n")

if (p_val_p < 0.05) {
  cat("Conclusion: The exponential rate differs significantly by prison history.\n")
} else {
  cat("Conclusion: No significant difference in time to exit by prison history.\n")
}
