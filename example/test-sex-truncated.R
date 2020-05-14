#### Simulating some data (liabilities) ####
h2 <- 0.5
K <- c(0.08, 0.02)
# K <- c(0.15, 0.05)
# K <- c(0.05, 0.05)
(thr <- qnorm(1 - K))

get_cov <- function(h2, n_sib = 0) {
  cov <- matrix(h2 / 2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3, 4] <- cov[4, 3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  cov
}

cov <- get_cov(h2, n_sib = 1)
liab <- MASS::mvrnorm(5000, rep(0, ncol(cov)), cov)

library(dplyr)
info_liab <- tibble(
  child_gen   = liab[, 1],
  child_full  = liab[, 2],
  child_sex   = sample(1:2, nrow(liab), replace = TRUE),
  mother_full = liab[, 3],
  mother_sex  = 1,
  father_full = liab[, 4],
  father_sex  = 2,
  sib_full    = liab[, 5],
  sib_sex     = sample(1:2, nrow(liab), replace = TRUE),
  #### Assign case-control status ####
  child_status  = (child_full  > thr[child_sex]),
  mother_status = (mother_full > thr[1]),
  father_status = (father_full > thr[2]),
  sib_status    = (child_full  > thr[sib_sex])
)

account_for_sex <- TRUE
if (!account_for_sex) thr[] <- qnorm(1 - mean(K))


#### LT-FH ####

nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, c(rep(list(c(FALSE, TRUE)), 3), list(1:2)))
names(all_config) <- c("child_status", "father_status", "mother_status", "child_sex")

# all_config <- do.call(expand.grid, c(rep(list(c(FALSE, TRUE)), nb_var), list(1:2, 1:2)))
# names(all_config) <- c("child_status", "father_status", "mother_status", "sib_status",
#                        "child_sex", "sib_sex")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config

set.seed(1)
simu_liab <- mvtnorm::rmvnorm(1e7, sigma = cov)

library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_sex     = sample(rep_len(1:2, nrow(simu_liab))),
  child_status  = (simu_liab[, 2] > thr[child_sex]),
  father_status = (simu_liab[, 3] > thr[2]),
  mother_status = (simu_liab[, 4] > thr[1]),
  sib_sex       = sample(rep_len(1:2, nrow(simu_liab))),
  sib_status    = (simu_liab[, 5] > thr[sib_sex])
) %>%
  left_join(all_config) %>%
  print()

group_means <- group_by(df_simu_liab, string, .drop = FALSE) %>%
  summarise(post_mean_liab = mean(child_gen), n = n(),
            se = sd(child_gen) / sqrt(n)) %>%
  ungroup() %>%
  arrange(post_mean_liab) %>%
  print(n = Inf)

#### Assign group posterior mean genetic liabilities to individuals ####
post_liab <- info_liab %>%
  # sample_frac(0.1, weight = ifelse(child_status, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(post_liab, child_status)

library(ggplot2)
ggplot(post_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, 
                 color = as.factor(child_sex)), alpha = 0.3) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(post_liab, c(cor(post_mean_liab, child_gen), cor(child_status, child_gen)))
# 1 sib: 0.4090318 0.3339000
# 1 sib / no sex: 0.4038593 0.3339000
# 0 sib: 0.3914486 0.3339000
# 0 sib / no sex: 0.3847673 0.3339000

#### Gibbs sampler for TMN distribution ####

source('example/gibbs-sampler.R')

post_liab$post_gen_liab <- NA
post_liab$post_gen_liab_se <- NA

for (i in 1:nrow(post_liab)) {
  print(i)
  lims <- cbind(upper = rep(Inf, ncol(cov)), lower = rep(-Inf, ncol(cov)))
  x <- post_liab[i, ]
  lims[2, x$child_status  + 1] <- thr[x$child_sex]
  lims[3, x$mother_status + 1] <- thr[x$mother_sex]
  lims[4, x$father_status + 1] <- thr[x$father_sex]
  lims[5, x$sib_status    + 1] <- thr[x$sib_sex]
  # fixed <- c(FALSE, x$child_status, x$mother_status, x$father_status)
  fixed <- rep(FALSE, ncol(cov))
  gen_liabs <- rtmvnorm.gibbs(20e3, burn_in = 1000, sigma = cov, 
                              lower = lims[, "lower"], 
                              upper = lims[, "upper"],
                              fixed = fixed)
  # gen_liabs <- gen_liabs[seq(1, length(gen_liabs), by = 20)]  # thinning
  post_liab$post_gen_liab[i] <- mean(gen_liabs)
  post_liab$post_gen_liab_se[i] <- sd(gen_liabs) / sqrt(length(gen_liabs))
}
summary(post_liab$post_gen_liab_se)

ggplot(post_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_gen_liab, child_gen, 
                 color = as.factor(child_sex)), alpha = 0.3) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(post_liab, c(cor(post_gen_liab, child_gen), cor(child_status, child_gen)))
# 0.4087596 0.3339000

ggplot(post_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_gen_liab, post_mean_liab, 
                 color = as.factor(child_sex)), alpha = 0.8) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")
