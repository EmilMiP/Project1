#### Simulate multivariate liabilities ####
h2 <- 0.5
cov <- diag(c(h2, 1 - h2, 1, 1))
colnames(cov) <- rownames(cov) <- 
  c("child_gen", "child_env", "father_full", "mother_full")
cov[3:4, 1] <- cov[1, 3:4] <- h2 / 2  # verify this but seems okay
cov

set.seed(1)
nsim <- 1e7
simu_liab <- mvtnorm::rmvnorm(nsim, sigma = cov)

is_case <- function(liab, age) {
  age_of_onset <- pnorm(liab, lower.tail = FALSE) * 500
  age > age_of_onset
}

library(dplyr)
df_simu_liab_train <- tibble(
  child_gen     = simu_liab[, 1],
  child_age     = runif(nsim, 10, 60),
  father_age    = child_age + runif(nsim, 20, 35),
  mother_age    = child_age + runif(nsim, 20, 35),
  child_status  = is_case(simu_liab[, 1] + simu_liab[, 2], child_age),
  father_status = is_case(simu_liab[, 3],                  father_age),
  mother_status = is_case(simu_liab[, 4],                  mother_age)
) %>%
  print()

colMeans(df_simu_liab_train)


mylm <- lm(child_gen ~ child_status + father_status + mother_status,
           data = df_simu_liab_train)
summary(mylm) # R2 = 0.1338 with 4 params

mylm <- lm(child_gen ~ child_status * (father_status + mother_status),
           data = df_simu_liab_train)
summary(mylm) # R2 = 0.1346 with 6 params // RSE = 0.6581


df_simu_liab_train$pred_lm <- predict(mylm, df_simu_liab_train)
# plot(pred_lm ~ post_mean_liab, data = sample_n(df_simu_liab_train, 50e3)); abline(0, 1, col = "red")

library(ggplot2)
ggplot(sample_n(df_simu_liab_train, 50e3)) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(pred_lm, child_gen), alpha = 0.3) + 
  geom_abline(col = "red") + 
  theme(legend.position = "top")

# unique(df_simu_liab_train$post_mean_liab) #  values
unique(df_simu_liab_train$pred_lm)        # 8 values


library(xgboost)
y <- df_simu_liab_train[[1]]
X <- as.matrix(df_simu_liab_train[2:7])
bst.params <- list(
  max_depth  = 2,
  base_score = 0,
  verbose    = 2,
  nthread    = 4
)
xgb <- xgboost(X, y, nrounds = 20, params = bst.params)

xgboost::xgb.importance(model = xgb)

df_simu_liab_train$pred_xgb <- predict(xgb, X)
# plot(pred_xgb ~ post_mean_liab, data = sample_n(df_simu_liab_train, 50e3)); abline(0, 1, col = "red")
plot(pred_xgb ~ pred_lm, data = sample_n(df_simu_liab_train, 50e3)); abline(0, 1, col = "red")

ggplot(sample_n(df_simu_liab_train, 50e3)) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(pred_xgb, child_gen), alpha = 0.3) + 
  geom_abline(col = "red") + 
  theme(legend.position = "top")

with(df_simu_liab_train, cor(pred_xgb, child_gen)^2)

