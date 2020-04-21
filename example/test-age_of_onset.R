#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 1e4
prev <- 0.05
downsample_frac <- 1 #1 for no downsampling.

set.seed(1) #set seed for more reliable comparison
father_gen <- rnorm(N, sd = sqrt(h2))
mother_gen <- rnorm(N, sd = sqrt(h2))
child_gen  <- (father_gen + mother_gen) / sqrt(2)

father_full <- father_gen + rnorm(N, sd = sqrt(1 - h2))
mother_full <- mother_gen + rnorm(N, sd = sqrt(1 - h2))
child_full  <- child_gen  + rnorm(N, sd = sqrt(1 - h2))

#### Assign case-control status ####
#K <- prev * 10:1/10
K <- prev * 5:1/5
#K <- prev * c(5,1)/5
#K <- c(0.05, 0.01)
#K <- prev
(thr <- qnorm(1 - K))


father_status <- rowSums(outer(father_full, thr, `>`))
mother_status <- rowSums(outer(mother_full, thr, `>`))
child_status  <- rowSums(outer(child_full, thr, `>`))

#### Simulate multivariate liabilities ####
cov <- diag(c(h2, 1 - h2, 1, 1))
colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
cov[3:4, 1] <- cov[1, 3:4] <- h2 / sqrt(2)  # verify this but seems okay
cov
round(100 * cov(cbind(child_gen, child_env = child_full - child_gen, 
                      father_full, mother_full)), 2)

nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, rep(list(0:length(thr)), nb_var))
names(all_config) <- c("child_status", "father_status", "mother_status")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config

set.seed(1)
simu_liab <- mvtnorm::rmvnorm(1e7, sigma = cov)

library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_status  = (rowSums(outer(simu_liab[, 1] + simu_liab[, 2], thr, `>`))),
  father_status = (rowSums(outer(simu_liab[, 3], thr, `>`))),
  mother_status = (rowSums(outer(simu_liab[, 4], thr, `>`)))
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
child_group <- tibble(child_gen, child_status, father_status, mother_status) %>%
  sample_frac(downsample_frac, weight = ifelse(child_status, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(child_group, child_status)

trues = child_group %>% group_by(string) %>% summarise(true_mean_liab = mean(child_gen)) %>% left_join(group_means)

library(ggplot2)
ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen), alpha = 0.2) + 
  geom_point(data = trues, aes(post_mean_liab, true_mean_liab), col = "red") +
  geom_abline(col = "red")
c(cor(child_group$post_mean_liab, child_group$child_gen),
  cor((child_group$child_status > 0) + 0L, child_group$child_gen))
### 5% cases + N = 5k
# k = 0.05:           0.4266394 0.3332564
# k = c(0.05, 0.01):  0.4296884 0.3332564
# k = 5:1/100         0.4313874 0.3332564
# k = 10:1/200        0.4319736 0.3332564

### 5% cases + N = 10k
# k = 0.05:           0.4337232 0.3183950
# k = c(0.05, 0.01):  0.4372176 0.3183950
# k = 5:1/100         0.4391433 0.3183950
# k = 10:1/200        0.4390453 0.3183950

### 50% cases + N = 5k
# k = 0.05:           0.7152612 0.6662278
# k = c(0.05, 0.01):  0.7255952 0.6662278
# k = 5:1/100         0.7286037 0.6662278
# k = 10:1/200        0.7305081 0.6662278

### 50% cases + N = 10k
# k = 0.05:           0.7148283 0.6697706
# k = c(0.05, 0.01):  0.7231299 0.6697706
# k = 5:1/100         0.7262446 0.6697706
# k = 10:1/200        0.7278923 0.6697706