#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 5000

father_gen <- rnorm(N, sd = sqrt(h2))
mother_gen <- rnorm(N, sd = sqrt(h2))
child_gen  <- (father_gen + mother_gen) / sqrt(2)

father_full <- father_gen + rnorm(N, sd = sqrt(1 - h2))
mother_full <- mother_gen + rnorm(N, sd = sqrt(1 - h2))
child_full  <- child_gen  + rnorm(N, sd = sqrt(1 - h2))

plot(father_gen, child_gen, pch = 20)
plot(mother_gen, child_gen, pch = 20)
round(100 * cor(cbind(father_gen, mother_gen, child_gen, 
                      father_full, mother_full, child_full))^2, 2)

#### Assign case-control status ####
K <- 0.1
(thr <- qnorm(1 - K))

father_status <- (father_full > thr)
mother_status <- (mother_full > thr)
child_status  <- (child_full  > thr)

#### Simulate multivariate liabilities ####
cov <- diag(c(h2, 1 - h2, 1, 1))
colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
cov[3:4, 1] <- cov[1, 3:4] <- h2 / sqrt(2)  # verify this but seems okay
cov
round(100 * cov(cbind(child_gen, child_env = child_full - child_gen, 
                      father_full, mother_full)), 2)

nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), nb_var))
names(all_config) <- c("child_status", "father_status", "mother_status")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config

simu_liab <- mvtnorm::rmvnorm(1e7, sigma = cov)

library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_status  = ((simu_liab[, 1] + simu_liab[, 2]) > thr),
  father_status = (simu_liab[, 3] > thr),
  mother_status = (simu_liab[, 4] > thr)
) %>%
  left_join(all_config) %>%
  print()

group_means <- group_by(df_simu_liab, string, .drop = FALSE) %>%
  summarise(post_mean_liab = mean(child_gen), n = n(),
            se = sd(child_gen) / sqrt(n)) %>%
  ungroup() %>%
  arrange(post_mean_liab) %>%
  print()

#### Assign group posterior mean genetic liabilities to individuals ####
child_group <- tibble(child_status, father_status, mother_status) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()

table(child_group$post_mean_liab)
library(ggplot2)
ggplot() +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(child_group$post_mean_liab, child_gen), alpha = 0.2) + 
  geom_abline(col = "red")
c(cor(child_group$post_mean_liab, child_gen),
  cor(child_status, child_gen))
