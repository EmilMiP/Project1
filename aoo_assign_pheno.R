setwd("C:/Users/FIUN7293/CommandStationBeta/Project1")
source("aoo_sampler.R")

assign_pheno = function(input,
                        h2 = 0.5,
                        child_aoo_col = "child_status",
                        age.col = "age",
                        dad_aoo_col = "dad_status",
                        mom_aoo_col = "mom_status",
                        age_parent = "age_parent",
                        sib_aoo_col_prefix = "sib",
                        sib_aoo_col_suffix = "status",
                        new_col_name = "aoo_ltfh",
                        thres,
                        s){
  #constructs a vector with column names to use:
  cols = c("FID", "IID", child_status_col)
  #includes any siblings:
  if (s > 0) {
    for (i in 1:s) {
      cols = c(cols, 
               paste(sib_status_col_prefix, i, sib_status_col_suffix, sep = ""))
    }  
  }
  #add parents:
  cols = c(cols, dad_status_col, mom_status_col)
  
  #extracts relevant info from input data
  data = input %>%
    select_at(vars(cols))
  
  gen_lia_est = sampler(thres = thres, h2 = h2, s = s)
  
  est_config = apply(gen_lia_est[,1:(4 + 2*s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  data_str = apply(data[,3:(6 + 2*s)], MARGIN = 1, FUN = function(x) paste(x, collapse = ""))
  data[[new_col_name]] = NA
  for (config in 1:length(est_config)) {
    data[[new_col_name]][data_str %in% est_config[config]] = gen_lia_est[config, "means"][[1]]
  }
  return(data)
}

