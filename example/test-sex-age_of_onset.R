#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 5e4
prev <- c(0.08, 0.02) * 4  
account_for_sex <- FALSE
#account_for_age <- FALSE
downsample_frac <- .4 #1 for no downsampling.
NA_frac = 1/3

assign_status = function(input, thr, acc_for_sex = account_for_sex) {
  len_thr = length(thr)
  N_input = length(input)
  input_age = sample(1:len_thr, N_input, replace = TRUE)
  input_status = rowSums(outer(input, thr, `>`)) + len_thr
  input_ctrl = which(input_status == len_thr)
  input_status[input_ctrl] = input_age[input_ctrl]
  return(input_status)
}
assign_status_sex = function(input, thr, input_sex, acc_for_sex = account_for_sex) {
#  N_input = length(input)9½2 input_sex = sample(1:2, N_input, replace = TRUE)
  input_status = c()
  input_status[input_sex == 1] = assign_status(input[input_sex == 1], thr[,1], acc_for_sex = account_for_sex)
  input_status[input_sex == 2] = assign_status(input[input_sex == 2], thr[,2], acc_for_sex = account_for_sex)
  return(input_status)
}

#K <- outer(10:1/10, prev, "*")
#K <- outer(5:1/5, prev, "*")
#K <- outer(c(5, 2)/5, prev, "*")
K <- matrix(prev,ncol = 2)
(thr <- qnorm(1 - K))

set.seed(1)
father_gen <- rnorm(N, sd = sqrt(h2))
mother_gen <- rnorm(N, sd = sqrt(h2))
child_gen  <- (father_gen + mother_gen) / sqrt(2)

father_full <- father_gen + rnorm(N, sd = sqrt(1 - h2))
mother_full <- mother_gen + rnorm(N, sd = sqrt(1 - h2))
child_full  <- child_gen  + rnorm(N, sd = sqrt(1 - h2))

child_sex <- sample(1:2, N, replace = TRUE)


father_status <- assign_status(father_full, thr[,1], acc_for_sex = TRUE)
mother_status <- assign_status(mother_full, thr[,2], acc_for_sex = TRUE)


#did not work, "cannot allocate vector of size 372.g Gb"
#child_status  <- rowSums(outer(child_full, thr[,child_sex], `>`))
child_status = assign_status_sex(child_full, thr, child_sex, acc_for_sex = TRUE)

if (!account_for_sex) thr <- matrix(qnorm(1 - apply(K, MARGIN = 1, mean)), ncol = 2, nrow = nrow(K)) #if all thresholds are the same, then some groups will not get any assigned to them
#if (!account_for_age) {
#  thr <- thr[1,]
#}

#### Simulate multivariate liabilities ####
cov <- diag(c(h2, 1 - h2, 1, 1))
colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
cov[3:4, 1] <- cov[1, 3:4] <- h2 / sqrt(2)  # verify this but seems okay


nb_var <- ncol(cov) - 1
all_config <- do.call(expand.grid, c(rep(list(c(NA,1:(nrow(K)*2))), nb_var), list(1:2))) ### change
names(all_config) <- c("child_status", "father_status", "mother_status", "child_sex")
all_config$string <- as.factor(do.call(paste, all_config + 0L))
all_config


set.seed(1)
simu_liab <- mvtnorm::rmvnorm(2e7, sigma = cov)


for (i in 2:ncol(cov)) { #starting at 2 to not filter on os twice, i.e. environment and genetic filtering
  simu_liab[base::sample(1:nrow(simu_liab), nrow(simu_liab) * NA_frac),i] <- NA
}


library(dplyr)
df_simu_liab <- tibble(
  child_gen     = simu_liab[, 1],
  child_sex     = sample(1:2, nrow(simu_liab), replace = TRUE),
  child_status  = assign_status_sex(simu_liab[, 1] + simu_liab[, 2], thr, child_sex),
  father_status = assign_status(simu_liab[, 3], thr[,1]),
  mother_status = assign_status(simu_liab[, 4], thr[,2])
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

child_group <- tibble(child_gen, child_status, father_status, mother_status, child_sex) %>%
  sample_frac(downsample_frac, weight = ifelse(child_status > nrow(thr), 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
count(child_group, child_status, child_sex)
count(child_group, child_status, "GWAS" = (child_status > nrow(thr)) + 0L)

trues = child_group %>% group_by(string) %>% summarise(true_mean_liab = mean(child_gen), n_true = n(), true_se = sd(child_gen) / sqrt(n_true)) %>% left_join(group_means)

library(ggplot2)
ggplot(child_group) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, child_gen, 
                 color = as.factor(child_sex)), alpha = 0.3) + 
#  geom_point(data = trues, aes(post_mean_liab, true_mean_liab)) +
  geom_abline(col = "black") + 
  theme(legend.position = "top")
with(child_group, c(cor(post_mean_liab, child_gen), cor(child_status > nrow(K), child_gen)))

####/w sex + age, 50% cases
#prev <- c(0.08, 0.02) + N = 10k      0.7161538 0.6231014

####/w age, 50% cases
#prev <- c(0.08, 0.02) + N = 10k      0.7271203 0.6362876


#### LTFH
# prev <- 0.05 , N = 10k,  50% cases  0.7098228 0.6697183

#### LTFH 
# prev <- c(0.05, 0.05) * 2, N = 10k, 50% cases  0.6985646 0.6340386

#all, 
#0.7195118 0.6530479 vs 0.6923003 0.6438172



