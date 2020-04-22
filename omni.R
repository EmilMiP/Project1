library(dplyr)

# Assisting Functions: ----------------------------------------------------
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





#### Simulating some data (liabilities) ####
h2 <- 0.5
N <- 5e4
prev <- c(0.08, 0.02) * 4  
account_for_sex <- FALSE
#account_for_age <- FALSE
downsample_frac <- .4 #1 for no downsampling.
NA_frac = 1/3
nsib <- 2

#### parameters needed:
#threshold, should reflect sex difference and aoo
#input data
###column names needed, i.e. sex column, name of parents etc.
#number of siblings to use
#fraction of NA values in simulated data
#number of simulated data in initial simulation
#heritability
#

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

if ( nsib > 0) { ### This handling of the siblings might be a major assumption, i.e. sharing genetic liability with the siblings
  sib_full <- as_tibble(
    sapply(1:nsib, FUN = function(x) child_gen + rnorm(N, sd = sqrt(1 - h2)))
  )
  colnames(sib_full) = paste("sib", 1:nsib ,"_full", sep = "")
  
  sib_sex <- as_tibble(
    sapply(1:nsib, FUN = function(x) sample(1:2, N, replace = TRUE))
  )
  colnames(sib_sex) <- paste("sib", 1:nsib, "_sex", sep = "")
}

input <- tibble(
  FID           = 1:N,
  IID           = 1:N,
  child_gen     = child_gen,
  child_sex     = child_sex,
  child_status  = assign_status_sex(child_full, thr, child_sex, acc_for_sex = TRUE),
  father_status = assign_status(father_full, thr[,1], acc_for_sex = TRUE),
  mother_status = assign_status(mother_full, thr[,2], acc_for_sex = TRUE),
)

if ( nsib > 0) {
  for (i in 1:nsib) {
    input[[paste("sib", i, "_sex", sep = "")]] <- sib_sex[[i]]
    input[[paste("sib", i, "_status", sep = "")]] <- assign_status_sex(sib_full[[i]], thr, sib_sex[[i]])
  }
}







omniscience = function(input, h2, thr, nsib = 0, NA_frac = .3, init_samp = 2e7, account_for_sex = FALSE, downsample_frac = 1) {

  ###assuming names are fixed for now, but might need to make them variable at a later point  
  
#### Simulate multivariate liabilities ####
if (nsib <= 0) {
  cov <- diag(c(h2, 1 - h2, 1, 1))
  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", "father_full", "mother_full")
  cov[2 + 1:2, 1] <- cov[1, 2 + 1:2] <- h2 / sqrt(2)  # verify this but seems okay
#  cov
#  round(100 * cov(cbind(child_gen, child_env = child_full - child_gen, 
#                        father_full, mother_full)), 2)
} else {### This handling of the siblings might be a major assumption
  cov <- diag(c(h2, 1 - h2, rep(1 - h2, nsib), 1, 1))
  colnames(cov) <- rownames(cov) <- c("child_gen", "child_env", if ( nsib > 0) paste("sib",1:nsib, "_full", sep = ""),"father_full", "mother_full")
  #fixing cov for gen child:
  cov[2 + nsib + 1:2, 1] <- cov[1, 2 + nsib + 1:2] <- h2 / sqrt(2)  # verify this but seems okay
  #fixing gen cov for siblings full and gen child:
  cov[2 + 1:nsib, 1] <- cov[1, 2 + 1:nsib] <- h2  # verify this
  # sorting out the correlations of the siblings
  cov[2 + 1:nsib, 2 + 1:nsib] <- cov[2 + 1:nsib, 2 + 1:nsib] + h2 # verify this 
  #sorting out the correlation between siblings and parents:
  cov[ncol(cov) - 0:1, ncol(cov) - 1 - 1:nsib] <- cov[ncol(cov) - 1 - 1:nsib, ncol(cov) - 0:1] <- h2 / sqrt(2) # verify this 
#  print(cov)
#  round(100 * cov(cbind(child_gen, 
#                        child_env = child_full - child_gen, 
#                        sib_full = sib_full,
#                        father_full, mother_full)), 2)
}

#this step becomes slow if the number of groups is high and sex is being accounted for (!)
nb_var <- ncol(cov) - 1
if (!account_for_sex) {
  all_config <- do.call(expand.grid, c(rep(list(c(NA,1:(nrow(thr)*2))), nb_var)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""), "father_status", "mother_status")
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
} else {
  all_config <- do.call(expand.grid, c(rep(list(c(NA,1:(nrow(thr)*2))), nb_var), rep(list(1:2), 1 + nsib)))
  names(all_config) <- c("child_status", if (nsib > 0) paste("sib",1:nsib, "_status", sep = ""), "father_status", "mother_status", "child_sex", if (nsib > 0) paste("sib",1:nsib, "_sex", sep = ""))
  all_config$string <- as.factor(do.call(paste, all_config + 0L))
}
all_config

set.seed(1)
simu_liab <- mvtnorm::rmvnorm(init_samp, sigma = cov)


for (i in 2:ncol(cov)) { #starting at 2 to not filter on os twice, i.e. environment and genetic filtering
  simu_liab[base::sample(1:nrow(simu_liab), nrow(simu_liab) * NA_frac),i] <- NA
}



if (nsib <= 0) {
  df_simu_liab <- tibble(
    child_gen     = simu_liab[, 1],
    child_sex     = sample(1:2, nrow(simu_liab), replace = TRUE),
    child_status  = assign_status_sex(simu_liab[, 1] + simu_liab[, 2], thr, child_sex),
    father_status = assign_status(simu_liab[, 3], thr[,1]),
    mother_status = assign_status(simu_liab[, 4], thr[,2])
  ) %>%
    left_join(all_config) %>%
    print()
} else {
  df_simu_liab <- tibble(
    child_gen     = simu_liab[, 1],
    child_sex     = sample(1:2, nrow(simu_liab), replace = TRUE),
    child_status  = assign_status_sex(simu_liab[, 1] + simu_liab[, 2], thr, child_sex),
    father_status = assign_status(simu_liab[, ncol(cov) - 1], thr[,1]),
    mother_status = assign_status(simu_liab[, ncol(cov)]    , thr[,2])
  ) 
  sibling_sex       = sapply(1:nsib, FUN = function(x) sample(1:2, nrow(simu_liab), replace = TRUE))
  sibling_status    = sapply(1:nsib, FUN = function(x) assign_status_sex(simu_liab[,2 + x], thr, sibling_sex[,x]))
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

child_group <- input %>%
  sample_frac(downsample_frac, weight = ifelse(input$child_status, 1e6, 1)) %>%
  left_join(all_config) %>%
  left_join(group_means) %>%
  print()
return(child_group)
}
