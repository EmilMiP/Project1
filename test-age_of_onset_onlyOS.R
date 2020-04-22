### implements the variable threshold for both cases and controls, rather than cases OR(!) controls.

#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 1e5
prev <- 0.2#c(0.08, 0.02)
downsample_frac <- .4 #1 for no downsampling.

assign_status = function(input, thr) { ## this function implicitly assumes the number of groups in case and controls are the same.
  len_thr = length(thr)
  N_input = length(input)
  input_age = sample(1:len_thr, N_input, replace = TRUE) ## age group is uniformly distributed
  input_status = rowSums(outer(input, thr, `>`)) + len_thr
  input_ctrl = which(input_status == len_thr)
  input_status[input_ctrl] = input_age[input_ctrl]
  return(input_status)
}

set.seed(1) #set seed for more reliable comparison
father_gen <- rnorm(N, sd = sqrt(h2))
mother_gen <- rnorm(N, sd = sqrt(h2))
child_gen  <- (father_gen + mother_gen) / sqrt(2)


child_full  <- child_gen  + rnorm(N, sd = sqrt(1 - h2))

#### Assign case-control status ####
#K <- prev * 10:1/10
K <- prev * 5:1/5
#K <- prev * c(5,1)/5
#K <- prev
(thr <- qnorm(1 - K))

len_thr = length(thr)

child_age = sample(1:len_thr, N, replace = TRUE)

child_status  <- assign_status(child_full, thr)


#### Simulate multivariate liabilities ####
cov <- diag(c(h2, 1 - h2))
colnames(cov) <- rownames(cov) <- c("child_gen", "child_env")
cov
round(100 * cov(cbind(child_gen, child_env = child_full - child_gen)), 2)

nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, rep(list(1:(len_thr*2)), nb_var))
names(all_config) <- c("child_status")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config



set.seed(1)
simu_liab <- mvtnorm::rmvnorm(1e6, sigma = cov)

library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_status  = assign_status(simu_liab[, 1] + simu_liab[, 2], thr),
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
child_group <- tibble(child_gen, child_status) %>%
  sample_frac(downsample_frac, weight = ifelse(child_status > len_thr, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(child_group, child_status)



trues = child_group %>% group_by(string) %>% summarise(true_mean_liab = mean(child_gen)) %>% left_join(group_means)

library(ggplot2)
ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen), alpha = 0.2) + 
  #  geom_point(data = trues, aes(post_mean_liab, true_mean_liab), col = "red") +
  geom_abline(col = "red")
c(cor(child_group$post_mean_liab, child_group$child_gen),
  cor((child_group$child_status > len_thr) + 0L, child_group$child_gen))


#0.4371630 0.3344324

#N = 5k, 50% cases
#0.7271625 0.6685495 vs 0.7137083 0.6738969

# N = 10k, 5% cases
#0.4525840 0.3445983 vs 0.4493351 0.3445983

# N = 100k, 50% cases
#0.7212471 0.6652323 vs  0.7113056 0.6652323

# N = 10k, 50% cases, prev = .1, 
# 0.7072824 0.6326771 vs 0.6946029 0.6326771, k = prev * c(5,1)/5
# 0.7263242 0.6453174 vs 0.6946029 0.6326771, k = prev * 5:1/5
# 0.7215895 0.6356807 vs 0.6946029 0.6326771, k = prev * 10:1/10
#I dont get why the correlation of the binary trait changes at all. can the cases be assigned differently each time ?

# N = 10k, 50% cases, prev = .2
# 0.7172276 0.5940128 vs 0.6911626 0.5940128, k = prev * c(5,1)/5
# 0.7302506 0.5985641 vs 0.6911626 0.5940128, k = prev * 5:1/5  ~ 5% improvement
# 0.7313964 0.5998677 vs 0.6911626 0.5940128, k = prev * 10:1/10