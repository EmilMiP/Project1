library(dplyr)
#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 1e5
account_for_sex <- T
K <- c(0.08, 0.02)
downsample_frac <- .1 #1 for no downsampling.
nsib <- 0
# K <- c(0.15, 0.05)

set.seed(1)
father_gen <- rnorm(N, sd = sqrt(h2))
mother_gen <- rnorm(N, sd = sqrt(h2))
child_gen  <- (father_gen + mother_gen) / sqrt(2)

child_full  <- child_gen  + rnorm(N, sd = sqrt(1 - h2))

child_sex <- sample(1:2, N, replace = TRUE)

#### Assign case-control status ####
(thr <- qnorm(1 - K))

child_status  <- (child_full  > thr[child_sex])
if (!account_for_sex) thr[] <- qnorm(1 - mean(K))


#### Simulate multivariate liabilities ####
cov <- diag(c(h2, 1 - h2))
colnames(cov) <- rownames(cov) <- c("child_gen", "child_env")
#cov[2 + 1:2, 1] <- cov[1, 2 + 1:2] <- h2 / sqrt(2)  # verify this but seems okay
cov
round(100 * cov(cbind(child_gen, child_env = child_full - child_gen)), 2)

nb_var <- ncol(cov) - 1
if (!account_for_sex) {
  all_config <- do.call(expand.grid, c(rep(list(c(FALSE, TRUE)), nb_var)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""))
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
} else {
  all_config <- do.call(expand.grid, c(rep(list(c(FALSE, TRUE)), nb_var), rep(list(1:2), 1 + nsib)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""), "child_sex", if (nsib > 0) paste("sib",1:nsib, "_sex", sep = ""))
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
}
all_config

set.seed(1)
simu_liab <- mvtnorm::rmvnorm(1e6, sigma = cov)


if (nsib <= 0) {
  df_simu_liab <- tibble(
    child_gen     = simu_liab[, 1],
    child_sex     = sample(1:2, nrow(simu_liab), replace = TRUE),
    child_status  = ((simu_liab[, 1] + simu_liab[, 2]) > thr[child_sex]),
  ) %>%
    left_join(all_config) %>%
    print()
} else {
  df_simu_liab <- tibble(
    child_gen     = simu_liab[, 1],
    child_sex     = sample(1:2, nrow(simu_liab), replace = TRUE),
    child_status  = ((simu_liab[, 1] + simu_liab[, 2]) > thr[child_sex]),
    father_status = (simu_liab[, ncol(cov) - 1] > thr[2]),
    mother_status = (simu_liab[, ncol(cov)] > thr[1])
  )
  sibling_sex       = sapply(1:nsib, FUN = function(x) sample(1:2, nrow(simu_liab), replace = TRUE))
  sibling_status    = sapply(1:nsib, FUN = function(x) simu_liab[,2 + x] > thr[sibling_sex[,x]])
  for (i in 1:nsib) {
    df_simu_liab[[paste("sib",i, "_sex", sep = "")]] <- sibling_sex[,i]
    df_simu_liab[[paste("sib",i, "_status", sep = "")]] <- sibling_status[,i]
  } #very inelegant, but it works..
  df_simu_liab = left_join(df_simu_liab, all_config) %>%
    print()
}

group_means <- group_by(df_simu_liab, string, .drop = FALSE) %>%
  summarise(post_mean_liab = mean(child_gen), n = n(),
            se = sd(child_gen) / sqrt(n)) %>%
  ungroup() %>%
  arrange(post_mean_liab) %>%
  print()

#### Assign group posterior mean genetic liabilities to individuals ####
if (nsib <= 0) {
  child_group <- tibble(child_gen, child_status, child_sex) %>%
    sample_frac(downsample_frac, weight = ifelse(child_status, 1e6, 1)) %>%
    left_join(all_config) %>%
    left_join(group_means) %>%
    print()
} else {
  child_group <- tibble(child_gen, child_status, father_status, mother_status, child_sex) 
  for (i in 1:nsib) {
    child_group[[paste("sib",i, "_sex", sep = "")]]    <- sib_sex[[i]]
    child_group[[paste("sib",i, "_status", sep = "")]] <- sib_status[,i]
  }
  
  child_group <- sample_frac(child_group, downsample_frac, weight = ifelse(child_status, 1e6, 1)) %>%
    left_join(all_config) %>%
    left_join(group_means) %>%
    print()
}

trues = child_group %>% group_by(string) %>% summarise(true_mean_liab = mean(child_gen)) %>% left_join(group_means)
library(ggplot2)
ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, 
                 color = as.factor(child_sex)), alpha = 0.3) +
  #  geom_point(data = trues, aes(post_mean_liab, true_mean_liab)) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(child_group, c(cor(post_mean_liab, child_gen), cor(child_status, child_gen)))
# K=c(0.08, 0.02) -> 0.4263897 0.3067237 vs 0.4214552 0.3067237
# K=c(0.15, 0.05) -> 0.5238736 0.3874479 vs 0.5171377 0.3874479
# with oversampling -> 0.6947441 0.6403549 vs 0.6916211 0.6415945
# K=c(0.08, 0.02) -> 0.4263897 0.3067237 vs 0.4214552 0.3067237
# K=c(0.15, 0.05) -> 0.5238736 0.3874479 vs 0.5171377 0.3874479
# with oversampling -> 0.6947441 0.6403549 vs 0.6916211 0.6415945
#0.7109676 0.6698215
