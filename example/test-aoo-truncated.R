#### Simulate multivariate liabilities ####
h2 <- 0.5
get_cov <- function(h2, n_sib = 0) {
  cov <- matrix(h2 / 2, 4 + n_sib, 4 + n_sib)
  diag(cov) <- 1
  cov[3, 4] <- cov[4, 3] <- 0
  cov[1:2, 1] <- cov[1, 1:2] <- h2
  cov
}

(cov <- get_cov(h2, n_sib = 2))

set.seed(1)
nsim <- 5000
simu_liab <- mvtnorm::rmvnorm(nsim, sigma = cov)

liab_to_aoo <- function(liab) pnorm(liab, lower.tail = FALSE) * 500

is_case <- function(liab, age) age > liab_to_aoo(liab)

aao_to_liab <- function(age) qnorm(age / 500, lower.tail = FALSE)
aao_to_liab(liab_to_aoo(0.15))

library(dplyr)
simu_liab <- tibble(
  child_gen   = simu_liab[, 1],
  child_age   = runif(nsim, 10, 60),
  father_age  = child_age + runif(nsim, 20, 35),
  mother_age  = child_age + runif(nsim, 20, 35),
  child_liab  = simu_liab[, 2],
  father_liab = simu_liab[, 3],
  mother_liab = simu_liab[, 4]
) %>%
  print()

colMeans(simu_liab)


#### Gibbs sampler for TMN distribution ####

source('example/gibbs-sampler.R')

simu_liab$post_gen_liab <- NA
simu_liab$post_gen_liab_se <- NA

for (i in 1:nrow(simu_liab)) {
  lower <- rep(-Inf, ncol(cov))
  upper <- rep(Inf, ncol(cov))
  x <- simu_liab[i, ]
  if (is_case(x$child_liab, x$child_age)) {
    lower[2] <- upper[2] <- x$child_liab
  } else {
    upper[2] <- aao_to_liab(x$child_age)
  }
  if (is_case(x$mother_liab, x$mother_age)) {
    lower[3] <- upper[3] <- x$mother_liab
  } else {
    upper[3] <- aao_to_liab(x$mother_age)
  }
  if (is_case(x$father_liab, x$father_age)) {
    lower[4] <- upper[4] <- x$father_liab
  } else {
    upper[4] <- aao_to_liab(x$father_age)
  }
  fixed <- (upper - lower) < 1e-4
  gen_liabs <- rtmvnorm.gibbs(20e3, burn_in = 1000, sigma = cov, 
                              lower = lower, upper = upper, fixed = fixed)
  # gen_liabs <- gen_liabs[seq(1, length(gen_liabs), by = 20)]  # thinning
  simu_liab$post_gen_liab[i] <- mean(gen_liabs)
  simu_liab$post_gen_liab_se[i] <- sd(gen_liabs) / sqrt(length(gen_liabs))
  batchmeans::bm(gen_liabs)
}
summary(simu_liab$post_gen_liab_se)

library(ggplot2)
ggplot(simu_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_gen_liab, child_gen, 
                 color = is_case(child_liab, child_age)), alpha = 0.3) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(simu_liab, c(cor(post_gen_liab, child_gen)**2, 
                  cor(is_case(child_liab, child_age)**2, child_gen)))
# 0.4497893 0.3538239

ggplot(simu_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_gen_liab, post_mean_liab, 
                 color = as.factor(child_sex)), alpha = 0.8) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")

#### LT-FH ####

nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), 3))
names(all_config) <- c("child_status", "father_status", "mother_status")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config

set.seed(1)
simu_ltfh <- mvtnorm::rmvnorm(1e7, sigma = cov)

library(dplyr)
df_simu_ltfh <- tibble(
  child_gen   = simu_ltfh[, 1],
  child_age   = runif(nrow(simu_ltfh), 10, 60),
  father_age  = child_age + runif(nrow(simu_ltfh), 20, 35),
  mother_age  = child_age + runif(nrow(simu_ltfh), 20, 35),
  child_status  = is_case(simu_ltfh[, 2], child_age),
  father_status = is_case(simu_ltfh[, 3], father_age),
  mother_status = is_case(simu_ltfh[, 4], mother_age)
) %>%
  left_join(all_config) %>%
  print()
mean(df_simu_ltfh$mother_status)
mean(df_simu_ltfh$child_status)

group_means <- group_by(df_simu_ltfh, string, .drop = FALSE) %>%
  summarise(post_mean_liab = mean(child_gen), n = n(),
            se = sd(child_gen) / sqrt(n)) %>%
  ungroup() %>%
  arrange(post_mean_liab) %>%
  print(n = Inf)

#### Assign group posterior mean genetic liabilities to individuals ####
post_liab <- simu_liab %>%
  mutate(child_status  = is_case(child_liab,  child_age),
         father_status = is_case(father_liab, father_age),
         mother_status = is_case(mother_liab, mother_age)) %>%
  # sample_frac(0.1, weight = ifelse(child_status, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(post_liab, child_status)

library(ggplot2)
ggplot(post_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, color = child_status), alpha = 0.3) + 
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(post_liab, c(cor(post_mean_liab, child_gen)**2, cor(child_status, child_gen)**2))
# 0.4326187 0.3538239

ggplot(post_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_gen_liab, child_gen, 
                 color = child_status), alpha = 0.3) + 
  geom_abline(col = "black") + 
  labs(x='Posterior Mean Liability (LTpred)', y='True genetic liability',color='Status')+
  scale_color_hue(labels = c("Control", "Case"))+
  theme(legend.position = "right")
with(post_liab, c(cor(post_gen_liab, child_gen)**2, cor(child_status, child_gen)**2))
# 0.4087596 0.3339000

ggplot(post_liab) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, 
                 color = child_status), alpha = 0.3) + 
  geom_abline(col = "black") + 
  labs(x='Posterior Mean Liability (LT-FH)', y='True genetic liability',color='Status')+
  scale_color_hue(labels = c("Control", "Case"))+
  theme(legend.position = "right")
with(post_liab, c(cor(post_mean_liab, child_gen)**2, cor(child_status, child_gen)**2))
# 0.4087596 0.3339000

bigstatsr::plot_grid(
  ggplot(post_liab) +
    bigstatsr::theme_bigstatsr() + 
    geom_point(aes(post_mean_liab, child_gen, 
                   color = child_status), alpha = 0.3) + 
    geom_abline(col = "black") + 
    labs(x='Posterior Mean Genetic Liability (LT-FH)', y='True genetic liability',color='Status')+
    scale_color_hue(labels = c("Control", "Case"))+
    theme(legend.position = "right")+
    xlim(-0.5, 2),
  ggplot(post_liab) +
    bigstatsr::theme_bigstatsr() + 
    geom_point(aes(post_gen_liab, child_gen, 
                   color = child_status), alpha = 0.3) + 
    geom_abline(col = "black") + 
    labs(x='Posterior Mean Genetic Liability (LTpred)', y='True genetic liability',color='Status')+
    scale_color_hue(labels = c("Control", "Case"))+
    theme(legend.position = "right")+
    xlim(-0.5, 2),
  
  ncol = 1
)

