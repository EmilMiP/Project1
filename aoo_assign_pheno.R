#setwd("C:/Users/FIUN7293/CommandStationBeta/Project1")
source("aoo_sampler.R")


input = phen
h2 = 0.5
child_aoo_col = "aoo"
age_col = "age"
dad_aoo_col = "dad_aoo"
mom_aoo_col = "mom_aoo"
age_parent = "age_parents"
sib_aoo_col_prefix = "sib"
sib_aoo_col_suffix = "_aoo"
new_col_name = "aoo_ltfh"
nthreads = 10
prev = .05
thres = aoo.thres
s = 2

assign_pheno = function(input,
                        h2 = 0.5,
                        child_aoo_col = "aoo",
                        age_col = "age",
                        dad_aoo_col = "dad_aoo",
                        mom_aoo_col = "mom_aoo",
                        age_parent = "age_parents",
                        sib_aoo_col_prefix = "sib",
                        sib_aoo_col_suffix = "_aoo",
                        new_col_name = "aoo_ltfh",
                        nthreads = 10,
                        prev = .05,
                        thres,
                        s){
  #constructs a vector with column names to use:
  cols = c("FID", "IID", age_col,  child_aoo_col)
  #includes any siblings:
  if (s > 0) {
    for (i in 1:s) {
      cols = c(cols, 
               paste(sib_aoo_col_prefix, i, sib_aoo_col_suffix, sep = ""))
    }  
  }
  #add parents:
  cols = c(cols, age_parent, dad_aoo_col, mom_aoo_col)
  
  #extracts relevant info from input data
  data = input %>%
    select_at(vars(cols))
  #assigning phenotype to people with aoo different than 0:
  data[[new_col_name]] = NA
  aoo.people = data[[child_aoo_col]] != 0
  #adding an artificial observation so the emperical cdf does not return 1 for the largest obs, meaning no Inf obs are generated.
  cum_prev = ecdf(c(data[[child_aoo_col]][aoo.people], max(data[[child_aoo_col]][aoo.people]) + 0.01 )) 
  #assigning the phenotype for cases, cased on their age of onset:
  data[[new_col_name]][data[[child_aoo_col]] != 0] = qnorm(1 - prev * (1 - cum_prev(data[[child_aoo_col]][aoo.people])), sd = sqrt(h2))
  
  
  gen_lia_est = sampler(thres = thres, h2 = h2, s = s)
  
  est_config = apply(gen_lia_est[,1:(3 + s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  data_str = apply(data[,grep("aoo", colnames(data))], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  data[[new_col_name]] = NA
  for (config in 1:length(est_config)) {
    data[[new_col_name]][data_str %in% est_config[config]] = gen_lia_est[config, "means"][[1]]
  }
  return(data)
}

plot( true$offspringgeno_lia[data[[child_aoo_col]] != 0], data[[new_col_name]][data[[child_aoo_col]] != 0])
abline(b = 1, a = 0)

cor(true$offspringgeno_lia[data[[child_aoo_col]] != 0], data[[new_col_name]][data[[child_aoo_col]] != 0])
