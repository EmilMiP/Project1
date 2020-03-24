library(truncnorm)
library(mvtnorm)
library(dplyr)
library(purrr)
library(stringr)
library(data.table)
library(doSNOW)
library(foreach)
library(parallel)
library(progress)

create.checker.mat = function(thres, s, no.age.grps = NULL, no.status = NULL){
  #"offspring_age" = 1:no.age.grps,
  expand.list = list("child_age" =  1:no.age.grps)
  expand.list[["child_grp"]] = 1:(no.status - 1) - 1
  if (s > 0) {
    for (i in 1:s) {
      expand.list[[paste("sib", i, "_age", sep = "")]] = 1:no.age.grps
      expand.list[[paste("sib", i, "_grp", sep = "")]] = 1:(no.status - 1) - 1
    }
  }
  #  expand.list[["parent_age"]] = 2:no.age.grps # start set to 2, since we dont want parents and offsping with the same age(except for the last grp.)
  expand.list[["dad_age"]] = 1:no.age.grps 
  expand.list[["dad_grp"]] = 1:(no.status - 1) - 1
  expand.list[["mom_age"]] = 1:no.age.grps
  expand.list[["mom_grp"]] = 1:(no.status - 1) - 1
  
  res = expand.grid(rev(expand.list))
  return(rev(res)) #rev used here to match the ordering of the grps with the output of dplyr
}

assign_obs = function(input, output, thresholds, s, checker.mat, output.col = "obs") {
  no.configs = dim(checker.mat)[1]
  grp.cols = grep("grp" ,colnames(checker.mat))
  age.cols = grep("age" ,colnames(checker.mat))
  
  input.dim = dim(input)
  rel.obs.ncol = input.dim[2] - 1
  rel.obs = matrix(NA, ncol = rel.obs.ncol, nrow = input.dim[1])
  for (i in 1:no.configs) {
    if (i %% 10 == 0) {cat("config:", i, "\n")}
    grps = checker.mat[i, grp.cols]
    ages = checker.mat[i, age.cols]
    ind = matrix(c(rep(ages, 2), grps + 1, grps + 2 ), ncol = 2)
    cur.thres = matrix(thresholds[ind], ncol = 2)
    
    for (ii in 1:rel.obs.ncol) {
      rel.obs[,ii] = input[, ii + 1] %between% cur.thres[ii ,]
    }
    
    output[[output.col]][[i]] =  c(output[[output.col]][[i]], input$child_raw[apply(rel.obs, MARGIN = 1, all)])
  }
  return(output)
}
sem <- function(x){
  if (length(x) > 1) {
    return(sd(x)/sqrt(length(x)))
  } else {
    return(Inf)
  }
}

resampler = function(sample_size = 5e5, cur.thres = cur.thres, s = s, h2 = h2, no.age.grps) {
  T_val_child = cur.thres[1,]
  T_val_dad = cur.thres[2 ,]
  T_val_mom = cur.thres[3,]
  
  #simulate parents lia. given status:
  dad_cond = rtruncnorm(n = sample_size,a = T_val_dad[1], b = T_val_dad[2] , mean = 0, sd = 1)
  mom_cond = rtruncnorm(n = sample_size,a = T_val_mom[1], b = T_val_mom[2] , mean = 0, sd = 1)
  
  #mean of (e_g, e_s) given (e_p1,e_p2)
  conditional_means = 0.5*h2*(dad_cond + mom_cond)
  
  #Construction of covariance matrix of (e_g, e_s) given (e_p1,e_p2):
  cov_matrix = diag(s + 1)
  cov_matrix[1,1] = h2*(1 - 0.5*h2)
  if (s > 0) {
    S = (diag(1 - 0.5*(h2)^2 - 0.5*h2*(1 - h2),s,s)) + (0.5*h2*(1 - h2))
    buffer = matrix(rep(0.5*h2*(1 - h2),s),nrow = s,ncol = 1)
    buffer_T = t(buffer)
    
    cov_matrix[1,2:(s + 1)] = buffer_T
    cov_matrix[2:(s + 1),1] = buffer
    cov_matrix[2:(s + 1),2:(s + 1)] = S
  }
  #Construction of covariance matrix of (e_g,e_s, e_p1, e_p2) given e_o (full lia. of offspring):
  cov_matrix_parents_eg_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),s + 3,s + 3)) + 0.5*h2*(1 - 0.5*h2)
  buffer_eg_sibs = matrix(rep(0.5*h2*(1 - h2),s + 2),nrow = s + 2,ncol = 1)
  cov_matrix_parents_eg_sibs[1,1] = h2*(1 - h2)
  cov_matrix_parents_eg_sibs[1,2:(s + 3)] = t(buffer_eg_sibs)
  cov_matrix_parents_eg_sibs[2:(s + 3),1] = buffer_eg_sibs
  cov_matrix_parents_eg_sibs[(s + 2),(s + 3)] = cov_matrix_parents_eg_sibs[(s + 3),(s + 2)]  = -1*(0.5*h2)^2
  
  # simulation of (e_g, e_s) given (e_p1, e_p2) for all configurations of parental stati:
  if (s > 0) {
    offspring_genetic_sibling_eps  = rmvnorm(n = sample_size, sigma = cov_matrix) + cbind(conditional_means, matrix(rep(conditional_means, s), ncol = s))
    sib.col.names = paste("sib", 1:s, sep = "")
    colnames(offspring_genetic_sibling_eps)[2:(1 + s)] = sib.col.names
  } else {
    offspring_genetic_sibling_eps  = rnorm(n = sample_size, mean = conditional_means, cov_matrix)
  }
  #first set of offspring simulations is done here.
  
  #simulate offsprings lia. given status:
  offspring_epsilon_o = rtruncnorm(n = sample_size, a = T_val_child[1], b = T_val_child[2], mean = 0, sd = 1)
  
  #calculating mean and simulating from (e_g, e_s, e_p1, e_p2) given e_o (full lia).
  #it is done in both configurations:
  offspring_genetic_sibling_parents_eo_means = matrix(c(h2*offspring_epsilon_o, rep(0.5*h2*offspring_epsilon_o,s + 2)), ncol = 3 + s )
  offspring_genetic_sibling_parents_eo = rmvnorm(n = sample_size, sigma = cov_matrix_parents_eg_sibs) + offspring_genetic_sibling_parents_eo_means
  if (s > 0) {
    colnames(offspring_genetic_sibling_parents_eo) = c("child_raw", sib.col.names, "dad", "mom")
  }
  #second set of sims are done here. No need to simulate environment, since it is implicitly done through the full lia.
  
  if ( s > 0) {
    T_val_sib = matrix(NA, ncol = 2, nrow = s)
    for (ii in 1:s) {
      T_val_sib[ii,] = cur.thres[ 3 + ii,]
    }
    
    fam_cond_sibs_status = cond_sibs(s = s, h2 = h2, T_val_sib = T_val_sib, sample_size = sample_size)
    no.rows = 2*sample_size + dim(fam_cond_sibs_status)[1]
  } else {
    no.rows = 2*sample_size
  }
  
  
  #the holder of stuff:
  df = data.frame("child_raw" = numeric(no.rows))
  if (s > 0) {
    #create the columns needed for the siblings:
    for (i in 1:s) {
      df[[paste("sib", i, sep = "")]] = numeric(no.rows)
    }
    
  }
  df[["dad"]] = numeric(no.rows) 
  df[["mom"]] = numeric(no.rows)
  df[["child"]] = numeric(no.rows) 
  
  #### Storing simulated values ####
  
  #note the order is important, since we want the observations to match each other.
  #e.g. the sims with offspring being a case needs to match
  
  #genetic liability for offspring:
  df$child_raw = c(offspring_genetic_sibling_eps[,1],
                   offspring_genetic_sibling_parents_eo[,1],
                   fam_cond_sibs_status$child_raw)
  
  fam_cond_sibs =  unname(fam_cond_sibs_status[,grep("sib", colnames(fam_cond_sibs_status))])
  names(fam_cond_sibs) = NULL
  #Liability for siblings:
  df[,paste("sib", 1:s, sep = "")] = rbind(offspring_genetic_sibling_eps[,2:(s + 1)],
                                           offspring_genetic_sibling_parents_eo[,2:(s + 1)],
                                           fam_cond_sibs_status[,grep("sib", colnames(fam_cond_sibs_status))])
  
  #Full liability for offspring:
  df$child = c(offspring_genetic_sibling_eps[,1] + rnorm(n = sample_size, mean = 0, sd = sqrt(1 - h2)),
               offspring_epsilon_o,
               fam_cond_sibs_status$child_raw + rnorm(n = sample_size, mean = 0, sd = sqrt(1 - h2)))
  
  #Full liability for dad:
  df$dad = c(dad_cond,
             offspring_genetic_sibling_parents_eo[,(s + 2)],
             fam_cond_sibs_status$dad)
  #Full liability for mom:
  df$mom = c(mom_cond,
             offspring_genetic_sibling_parents_eo[,(s + 3)],
             fam_cond_sibs_status$mom)
  
  
  #df = assign_grps_resamp(input = df, thresholds = thres, s = s, no.age.grps = no.age.grps)
  #
  ## cols = grep(paste(c( "sex", "status"), collapse = "|"), colnames(df))
  #grp.means.resamp = df %>%
  #  group_by_at(vars(cols)) %>%
  #  summarise("obs" = list(child_raw)) %>%
  #  ungroup()
  return(df)
}

cond_sibs = function(s, h2, T_val_sib, sample_size) {
  rho = 0.5*h2
  
  if (s == 1) {
    sib = rtruncnorm(n = sample_size, a = T_val_sib[1,1], b = T_val_sib[1,2])
    cov_cond_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),3,3)) + 0.5*h2*(1 - 0.5*h2)
    cov_cond_sibs[1,1] = h2*(1 - h2)
    cov_cond_sibs[2, 3] = cov_cond_sibs[3,2]  = -1*(0.5*h2)^2
    fam_cond_sib = rmvnorm(n = sample_size, mean = rep(0,3), sigma = cov_cond_sibs) + matrix(rho*rep(sib, each = 3), byrow = T, ncol = 3)
    no.rows = sample_size
    df_fam = data.frame("child_raw" = numeric(no.rows),
                        "dad" = numeric(no.rows),
                        "mom" = numeric(no.rows),
                        "sib1" = numeric(no.rows))
    df_fam[,c("child_raw", "dad", "mom")] = fam_cond_sib
    df_fam[,"sib1"] = sib
    return(df_fam)
  }
  if (s == 2) {
    sib1 = rtruncnorm(n = sample_size, a = T_val_sib[1,1], b = T_val_sib[1,2])
    sib2 = rtruncnorm(n = sample_size, a = T_val_sib[2,1], b = T_val_sib[2,2])
    
    cov_cond_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),3,3)) + 0.5*h2*(1 - 0.5*h2)
    cov_cond_sibs[1, 1] = h2*(1 - h2)
    cov_cond_sibs[2, 3] = cov_cond_sibs[3,2]  = -1*(0.5*h2)^2
    
    fam_cond_sib = rmvnorm(n = sample_size, mean = rep(0,3), sigma = cov_cond_sibs) + matrix(rho*rep((sib1 + sib2), each = 3), byrow = T, ncol = 3)
    no.rows = sample_size
    df_fam = data.frame("child_raw" = numeric(no.rows),
                        "dad" = numeric(no.rows),
                        "mom" = numeric(no.rows),
                        "sib1" = numeric(no.rows),
                        "sib2" = numeric(no.rows))
    df_fam[,c("child_raw", "dad", "mom")] = fam_cond_sib
    df_fam[,"sib1"] = sib1
    df_fam[,"sib2"] = sib2
    return(df_fam)
  }
  if (s >= 3) {
    cov_cond_sibs = (diag(1 - (0.5*h2)^2 - 0.5*h2*(1 - 0.5*h2),3,3)) + 0.5*h2*(1 - 0.5*h2)
    cov_cond_sibs[1, 1] = h2*(1 - h2)
    cov_cond_sibs[2, 3] = cov_cond_sibs[3,2]  = -1*(0.5*h2)^2
    trunc_sibs = list()
    for (ii in 1:s) {
      trunc_sibs[[paste("sib", ii, sep = "")]] = rtruncnorm(n = sample_size, a = T_val_sib[ii,1], b = T_val_sib[ii,2])
    }
    trunc_sibs = as.data.frame(trunc_sibs)
    mu_trunc_sibs = rowSums(trunc_sibs) * rho * (1 - rho) / (1 + rho - 2*rho^2)
    fam_cond_sibs = rmvnorm(n = sample_size, mean = rep(0,3), sigma = cov_cond_sibs) + matrix(rep(mu_trunc_sibs, each = 3), byrow = T, ncol = 3)
    no.rows = sample_size
    df_fam = data.frame("child_raw" = numeric(no.rows),
                        "dad" = numeric(no.rows),
                        "mom" = numeric(no.rows))
    df_fam[,c("child_raw", "dad", "mom")] = fam_cond_sibs
    for (ii in 1:s) {
      df_fam[[paste("sib",ii, sep = "")]] = trunc_sibs[[ii]]
    }
    return(df_fam)
  }
}


h2 = .5
s = 2
k = .05
no.grps = 1
prev.thres = k*no.grps:1/no.grps
aoo.thres = matrix(NA, ncol = 3, nrow = length(prev.thres))

for (j in 1:(length(prev.thres))) {
  #aoo.thres[[1]][j] = list(c(-Inf, qnorm(1 - prev.thres[j])))
  aoo.thres[j,] = c(-Inf, qnorm(1 - prev.thres[j]), Inf)
} 
thres = aoo.thres

checker.mat = as.matrix(create.checker.mat(thres, s = 2, no.age.grps = dim(thres)[1], no.status = dim(thres)[2]))



no.age.grps = dim(thres)[1]
cov_matrix = diag(1 - 0.5*h2, s + 3, s + 3) + 0.5*h2 
cov_matrix[1,1] = h2
cov_matrix[(s + 2),(s + 3)] = cov_matrix[(s + 3),(s + 2)] = 0

child.env = rnorm(n = 5000000, sd = sqrt(1 - h2))
sim.dat = as.data.frame(rmvnorm(n = 5000000, mean = rep(0, 3 + s), sigma = cov_matrix))
if (s > 0) {
  colnames(sim.dat) = c("child_raw", paste("sib", 1:s, sep = ""), "dad", "mom")
} else {
  colnames(sim.dat) = c("child_raw", "dad", "mom")
}
#filling out values for the "child" column, i.e. full child liability.
sim.dat[["child"]] = sim.dat$child_raw + child.env 
#move columns to be in the same order as the grp means columns:
sim.dat[,2:(4 + s)] = sim.dat[,c("child", colnames(sim.dat)[-grep("child", colnames(sim.dat))])]
#moving names to reflect the columns:
colnames(sim.dat)[2:(4 + s)] = c("child", colnames(sim.dat)[-grep("child", colnames(sim.dat))])


#updated the checker mat for a new threshold format.
grp.means = as_tibble(create.checker.mat(thres, s = 2, no.age.grps = dim(thres)[1], no.status = dim(thres)[2]))
grp.means$obs <- list(NULL)

grp.means = assign_obs(input = sim.dat, output = grp.means, thresholds = thres, s, checker.mat = checker.mat)

sem.vec = sapply(grp.means$obs, FUN = sem)
if (any(sem.vec > 0.01)) {
  grp.means[["resamp_obs"]] = list(NULL)
}

while (any(sem.vec > 0.01)) {
  cat("configurations requiring resampling:" , sum(sem.vec > 0.01), "\n"  )
  ind.max.sem = which.max(sem.vec) ## maybe only take one at a time and update inbetween instead of all at once ?
  # entries = sample(which(sem.vec > 0.01), size = nthreads)
  # grp.means[entries,]
  # age.cols = grp.means[ind.max.sem, str_subset(colnames(grp.means), "age")]
  grps = checker.mat[ind.max.sem, grep("grp", colnames(checker.mat))]
  ages = checker.mat[ind.max.sem, grep("age", colnames(checker.mat))]
  ind = matrix(c(rep(ages, 2), grps + 1, grps + 2 ), ncol = 2)
  cur.thres = matrix(thres[ind], ncol = 2)
  # ph = foreach(i = 1:nthreads,
  #              .packages = c("data.table", "truncnorm", "mvtnorm", "dplyr")) %dopar% {
  grp.means.resamp = resampler(sample_size = 50000, cur.thres = cur.thres, s = s, h2 = h2, no.age.grps = no.age.grps)
  grp.means = assign_obs(input = grp.means.resamp, grp.means, thresholds = thres, s = s, checker.mat = checker.mat, output.col = "resamp_obs")
  #             }
  #  for (ii in seq_along(ph)) {
  #    update_entries = apply(ph[[ii]][,1:(3 + s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  #    #    match.ph = match.ph[!is.na(match.ph)]
  #    if (!all(grp.means.str[grp.means.str %in% update_entries] == update_entries)) {
  #      print("STOOOOOP")
  #    }
  #    update_entries_no = match(update_entries,grp.means.str)
  #    grp.means$obs[update_entries_no] = map2(grp.means$obs[update_entries_no], ph[[ii]]$obs, ~ c(.x, .y))
  #    
  #  }
  # update_entries =  apply(grp.means.resamp[,1:(3 + s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  # if (!all(grp.means.str[grp.means.str %in% update_entries] == update_entries)) {
  #   print("STOOOOOP")
  # }
  # update_entries_no = match(update_entries,grp.means.str)
  # grp.means$obs[update_entries_no] = map2(grp.means$obs[update_entries_no], grp.means.resamp$obs, ~ c(.x, .y))
  # #    grp.means[["mean"]] = map2_dbl(grp.means$all_obs, grp.means$all_resamp, ~ mean(c(.x,.y)))
  sem.vec = sapply(grp.means$obs, FUN = sem)
}
#stopCluster(cl)
grp.means[["means"]] = sapply(grp.means$obs, FUN = mean)
grp.means[["resamp_means"]] = sapply(grp.means$resamp_obs, FUN = mean)



par(mfrow = c(2,1))
hist(grp.means$obs[[32]], breaks = 50)
hist(grp.means$resamp_obs[[4]], breaks = 50)

head(sim.dat)
df = matrix(NA, ncol = 3 + s, nrow = dim(sim.dat)[1])
for (k in 1:(3 + s)) {
  df[,k] = sim.dat[,1 + k] %between% thres[-1]
}
#apply(df, MARGIN = 1, all))
mean(sim.dat[rowSums(df) == (3 + s),1])
