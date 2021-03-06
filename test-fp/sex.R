#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 1e4
account_for_sex <- TRUE
K <- c(0.08, 0.02)
# K <- c(0.15, 0.05)
# K <- c(0.05, 0.05)

set.seed(1)
father_gen <- rnorm(N, sd = sqrt(h2))
mother_gen <- rnorm(N, sd = sqrt(h2))
child_gen  <- (father_gen + mother_gen) / sqrt(2)

father_full <- father_gen + rnorm(N, sd = sqrt(1 - h2))
mother_full <- mother_gen + rnorm(N, sd = sqrt(1 - h2))
child_full  <- child_gen  + rnorm(N, sd = sqrt(1 - h2))

child_sex <- sample(1:2, N, replace = TRUE)

#### Assign case-control status ####
(thr <- qnorm(1 - K))

father_status <- (father_full > thr[2])
mother_status <- (mother_full > thr[1])
child_status  <- (child_full  > thr[child_sex])

if (!account_for_sex) thr[] <- qnorm(1 - mean(K))

#### Simulate multivariate liabilities ####
cov <- diag(c(h2, 1 - h2, 1, 1))
colnames(cov) <- rownames(cov) <- 
  c("child_gen", "child_env", "father_full", "mother_full")
cov[3:4, 1] <- cov[1, 3:4] <- h2 / 2  # verify this but seems okay
cov
round(100 * cov(cbind(child_gen, child_env = child_full - child_gen, 
                      father_full, mother_full)), 2)

nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, c(rep(list(c(FALSE, TRUE)), nb_var), list(1:2)))
names(all_config) <- c("child_status", "father_status", "mother_status", "child_sex")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config

set.seed(1)
simu_liab <- mvtnorm::rmvnorm(1e7, sigma = cov)

library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_sex     = rep_len(1:2, nrow(simu_liab)),
  child_status  = ((simu_liab[, 1] + simu_liab[, 2]) > thr[child_sex]),
  father_status = (simu_liab[, 3] > thr[2]),
  mother_status = (simu_liab[, 4] > thr[1])
) %>%
  mutate(child_sex_bool = (child_sex > 1.5)) %>%
  left_join(all_config) %>%
  print()

group_means <- group_by(df_simu_liab, string, .drop = FALSE) %>%
  summarise(post_mean_liab = mean(child_gen), n = n(),
            se = sd(child_gen) / sqrt(n)) %>%
  ungroup() %>%
  arrange(post_mean_liab) %>%
  print()

#### Assign group posterior mean genetic liabilities to individuals ####
child_group <- tibble(child_gen, child_status, father_status, mother_status, child_sex) %>%
  # sample_frac(0.1, weight = ifelse(child_status, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(child_group, child_status)

library(ggplot2)
ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, 
                 color = as.factor(child_sex)), alpha = 0.3) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(child_group, c(cor(post_mean_liab, child_gen), cor(child_status, child_gen)))
# K=c(0.08, 0.02) -> 0.4263897 0.3067237 vs 0.4214552 0.3067237
# K=c(0.15, 0.05) -> 0.5238736 0.3874479 vs 0.5171377 0.3874479
# with oversampling -> 0.6947441 0.6403549 vs 0.6916211 0.6415945

mylm <- lm(child_gen ~ child_status + father_status + mother_status,
           data = df_simu_liab)
summary(mylm) # 0.1419 with 4 params

mylm <- lm(child_gen ~ child_status + father_status + mother_status +
             child_status:father_status + child_status:mother_status,
           data = df_simu_liab)
summary(mylm) # 0.1425 with 6 params

mylm <- lm(child_gen ~ child_status * father_status * mother_status * child_sex,
           data = df_simu_liab)
summary(mylm) # 0.1451 with 16 params

mylm <- lm(child_gen ~ child_status * I(father_status + mother_status) * child_sex,
           data = df_simu_liab)
summary(mylm) # 0.1448 with 8 params

mylm <- lm(child_gen ~ child_status * (father_status + mother_status) * child_sex,
           data = df_simu_liab)
summary(mylm) # 0.1451 with 12 params

mylm <- lm(child_gen ~ child_status * (father_status + mother_status) * child_sex +
             child_status:father_status:mother_status,
           data = df_simu_liab)
summary(mylm) # 0.1451 with 14 params


child_group$pred_lm <- predict(mylm, child_group)
plot(pred_lm ~ post_mean_liab, data = child_group); abline(0, 1, col = "red")

ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(pred_lm, child_gen, 
                 color = as.factor(child_sex)), alpha = 0.3) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")

unique(child_group$post_mean_liab) # 16 values
unique(child_group$pred_lm)        # 16 values


library(xgboost)
y <- df_simu_liab[[1]]
X <- as.matrix(df_simu_liab[2:5])
xgb <- xgboost(X, y, nrounds = 20, verbose = 2, nthreads = 4)
xgboost::xgb.importance(model = xgb)
xgboost::xgb.plot.tree(modle = xgb)
print(xgb)

child_group$pred_xgb <- predict(xgb, as.matrix(child_group[c(5, 2:4)]))
plot(pred_xgb ~ post_mean_liab, data = child_group); abline(0, 1, col = "red")
plot(pred_xgb ~ pred_lm, data = child_group); abline(0, 1, col = "red")

unique(child_group$pred_xgb)        # 16 values
