library(dplyr)
library(mvtnorm)
#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 1e4
account_for_sex <- F
multiplier = 1
K <- c(0.1, 0.1) * multiplier
downsample_frac <- 1 #.1 * multiplier #1 for no downsampling.
nsib <- 0
# K <- c(0.15, 0.05)

#### Simulate multivariate liabilities ####
if (nsib <= 0) {
  cov <- diag(c(h2, 1 - h2, 1, 1))
  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
  cov[2 + 1:2, 1] <- cov[1, 2 + 1:2] <- h2 / 2
  cov
} else {### This handling of the siblings might be a major assumption
  cov <- diag(c(h2, 1 - h2, rep(1 - h2 / 2, nsib), 1, 1))
  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", if ( nsib > 0) paste("sib",1:nsib, "_full", sep = ""),"father_full", "mother_full")
  #fixing cov for gen child:
  cov[2 + nsib + 1:2, 1] <- cov[1, 2 + nsib + 1:2] <- h2 / 2
  #fixing gen cov for siblings full and gen child:
  cov[2 + 1:nsib, 1] <- cov[1, 2 + 1:nsib] <- h2 / 2
  # sorting out the correlations of the siblings
  cov[2 + 1:nsib, 2 + 1:nsib] <- cov[2 + 1:nsib, 2 + 1:nsib] + h2 / 2 
  #sorting out the correlation between siblings and parents:
  cov[ncol(cov) - 0:1, ncol(cov) - 1 - 1:nsib] <- cov[ncol(cov) - 1 - 1:nsib, ncol(cov) - 0:1] <- h2 / 2 
  cov
}


set.seed(1)
family = rmvnorm(n = N, mean = rep(0, 4 + nsib), sigma = cov)
colnames(family) <- colnames(cov)

input <- tibble(
  child_gen    = family[,1],
  child_full   = family[,1] + family[,2],
  child_sex    = sample(1:2, N, replace = TRUE),
  father_full  = family[,ncol(family) - 1],
  mother_full  = family[, ncol(family)]
)


if ( nsib > 0) { ### This handling of the siblings might be a major assumption, i.e. sharing genetic liability with the siblings
  for (i in 1:nsib) {
    input[[paste("sib", i ,"_full", sep = "")]] = family[,2 + i]
    input[[paste("sib", i, "_sex", sep = "")]] = sample(1:2, N, replace = TRUE)
  }
}

#### Assign case-control status ####
(thr <- qnorm(1 - K))

input[["father_status"]] <- (input$father_full > thr[2])
input[["mother_status"]] <- (input$mother_full > thr[1])
input[["child_status"]]  <- (input$child_full  > thr[input$child_sex])
if ( nsib > 0) {
  for (i in 1:nsib) {
    input[[paste("sib", i, "_status", sep = "")]] = input[[4 + 2 * i]] > thr[input[[5 + 2 * i]]]
  }
}
if (!account_for_sex) thr[] <- qnorm(1 - mean(K))




nb_var <- ncol(cov) - 1
if (!account_for_sex) {
  all_config <- do.call(expand.grid, c(rep(list(c(FALSE, TRUE)), nb_var)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""),"father_status", "mother_status")
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
} else {
  all_config <- do.call(expand.grid, c(rep(list(c(FALSE, TRUE)), nb_var), rep(list(1:2), 1 + nsib)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""), "father_status", "mother_status", "child_sex", if (nsib > 0) paste("sib",1:nsib, "_sex", sep = ""))
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
}
all_config

set.seed(1)
simu_liab <- mvtnorm::rmvnorm(1e7, sigma = cov)


if (nsib <= 0) {
  df_simu_liab <- tibble(
    child_gen     = simu_liab[, 1],
    child_sex     = sample(1:2, nrow(simu_liab), replace = TRUE),
    child_status  = ((simu_liab[, 1] + simu_liab[, 2]) > thr[child_sex]),
    father_status = (simu_liab[, ncol(cov) - 1] > thr[2]),
    mother_status = (simu_liab[, ncol(cov)]     > thr[1])
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
  child_group <- input %>%
    sample_frac(downsample_frac, weight = ifelse(child_status, 1e6, 1)) %>%
    left_join(all_config) %>%
    left_join(group_means) %>%
    print()
} else {
  child_group <- input %>%
    sample_frac(downsample_frac, weight = ifelse(child_status, 1e6, 1)) %>%
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
