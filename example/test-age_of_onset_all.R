### implements the variable threshold for both cases and controls, rather than cases OR(!) controls.

#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 5e4
prev <- 0.2 #c(0.08, 0.02)
downsample_frac <- .4 #1 for no downsampling.
nsib = 0


assign_status = function(input, thr) { ## this function implicitly assumes the number of groups in case and controls are the same.
  len_thr = length(thr)
#  N_input = length(input)
#  input_age = sample(1:len_thr, N_input, replace = TRUE) ## age group is uniformly distributed
  input_status = rowSums(outer(input, thr, `>`)) + len_thr
  input_ctrl = which(input_status == len_thr)
  #these individuals did not pass any threshold, assigning lowest grp. 0 reserved for NA.
  input_status[input_ctrl] = 1
  return(input_status)
}
library(mvtnorm)





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
  ID           = 1:N,
  child_gen    = family[,1],
  child_full   = family[,1] + family[,2],
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
#K <- prev * 10:1/10
K <- prev * 5:1/5
#K <- prev * c(5,1)/5
K <- prev
(thr <- qnorm(1 - K))

len_thr = length(thr)
#### Assign case-control status ####

input[["father_status"]] <- assign_status(input$father_full, thr)
input[["mother_status"]] <- assign_status(input$mother_full, thr)
input[["child_status"]]  <- assign_status(input$child_full, thr)
if ( nsib > 0) {
  for (i in 1:nsib) {
    input[[paste("sib", i, "_status", sep = "")]] = input[[4 + 2 * i]] > thr[input[[5 + 2 * i]]]
  }
}

#### Simulate multivariate liabilities ####
cov <- diag(c(h2, 1 - h2, 1, 1))
colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
cov[3:4, 1] <- cov[1, 3:4] <- h2 / 2#sqrt(2)  # verify this but seems okay
cov


nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, rep(list(1:(len_thr*2)), nb_var))
names(all_config) <- c("child_status", "father_status", "mother_status")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config



set.seed(1)
simu_liab <- mvtnorm::rmvnorm(3e6, sigma = cov)

library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_status  = assign_status(simu_liab[, 1] + simu_liab[, 2], thr),
  father_status = assign_status(simu_liab[, 3], thr),
  mother_status = assign_status(simu_liab[, 4], thr)
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
child_group <- input %>%
  sample_frac(downsample_frac, weight = ifelse(child_status > len_thr, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(child_group, child_status)




library(ggplot2)
ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen), alpha = 0.2) + 
#  geom_point(data = trues, aes(post_mean_liab, true_mean_liab), col = "red") +
  geom_abline(col = "red")



input <- tibble(
  ID           = 1:N,
  child_gen    = family[,1],
  child_full   = family[,1] + family[,2],
  child_sex    = sample(1:2, N, replace = TRUE),
  father_full  = family[,ncol(family) - 1],
  mother_full  = family[, ncol(family)]
)


if ( nsib > 0) { ### This handling of the siblings might be a major assumption, i.e. sharing genetic liability with the siblings
  for (i in 1:nsib) {
    input[[paste("sib", i ,"_full", sep = "")]] = family[,2 + i]
    input[[paste("sib", i, "_sex", sep = "")]]  = sample(1:2, N, replace = TRUE)
  }
}
#prev = .1

res = matrix(NA, ncol = 3, nrow = 4)

Klist = list(
  prev,
  prev * c(5,1)/5,
  prev * 5:1/5
#  prev * 10:1/10
  )

for (i in seq_along(Klist)) {
  (thr <- qnorm(1 - Klist[[i]]))
  
#  len_thr = length(thr)
  #### Assign case-control status ####
  input[["father_status"]] <- assign_status(input$father_full, thr)
  input[["mother_status"]] <- assign_status(input$mother_full, thr)
  input[["child_status"]]  <- assign_status(input$child_full, thr)
  if ( nsib > 0) {
    for (i in 1:nsib) {
      input[[paste("sib", i, "_status", sep = "")]] = input[[4 + 2 * i]] > thr[input[[5 + 2 * i]]]
    }
  }
  
  thr = cbind(thr,thr)
  
  omni_dat = omniscience(input = input, thr = thr, h2 = 0.5, init_samp = 2e6, account_for_sex = FALSE, NA_frac = 1/(1 + length(thr)*2), nsib = nsib)
  colnames(omni_dat)[11] = "omni"
  ph = child_group %>% left_join(omni_dat, by = c("ID","child_gen"))
  res[i,] = c(cor(ph$post_mean_liab, ph$child_gen),
             cor(ph$child_status.x, child_group$child_gen),
             cor(ph$omni, ph$child_gen))
}

res
res[,3]/res[,1]

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

#0.6724999 0.5300688 0.6764113
#0.6527389 0.5300688 0.6673333