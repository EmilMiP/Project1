options(scipen = 999)
library(data.table)
library(stringr)
library(dplyr)
source("D:\\Work\\Project1\\Code\\Project1\\omni.R")

#get list of all the summary stats to use, and their true values.
#fileRoot = "/home/emp/faststorage/project1/simulatedData/"
fileRoot = "D:\\Work\\Project1\\SimulatedData\\"
base = "sibs2"
sword = "100kx10k"
cword = "C1000"
bword = "NF"
identifier = "maf005"
#exword = 


true.files = list.files(path = fileRoot, pattern = "true", full.names = T)
true.files = str_subset(true.files, pattern = base)
true.files = str_subset(true.files, pattern = sword)
true.files = str_subset(true.files, pattern = cword)
true.files = str_subset(true.files, pattern = bword)
true.files = str_subset(true.files, pattern = identifier)


phen.files = list.files(path = fileRoot, pattern = "\\.phen", full.names = T)
phen.files = str_subset(phen.files, pattern = base)
phen.files = str_subset(phen.files, pattern = sword)
phen.files = str_subset(phen.files, pattern = cword)
phen.files = str_subset(phen.files, pattern = bword)
phen.files = str_subset(phen.files, pattern = identifier)
phen.files = str_subset(phen.files, pattern = "\\\\sibs")

# phen.files = phen.files[(0:9)*2+ 1]
multiplier = 4
prev = c(.08, 0.02) * multiplier


## updating phenotype file for ltfh to reflect the sex differences
for (i in seq_along(true.files)) {
  thr = qnorm(1 - prev)
  phen = fread(phen.files[i])
  true = fread(true.files[i])
  res = as.data.frame(matrix(NA, nrow = nrow(phen), ncol = ncol(phen)))
  res[,1] = phen[[1]]
  res[,2] = phen[[2]]
  colnames(res) = colnames(phen)
  lia_indiv_vec = colnames(true)[grep("_lia", colnames(true))]
  lia_indiv_vec = lia_indiv_vec[-grep("geno", lia_indiv_vec)]
  sex_indiv_vec = colnames(true)[grep("sex", colnames(true))]
  sex_indiv_vec[4:5] = sex_indiv_vec[2:3]
  sex_indiv_vec[2:3] = c("father", "mother")
  ctr = 1
  for (ii in 1:length(lia_indiv_vec)) {
    lia_indiv = lia_indiv_vec[ii]
    sex_indiv = sex_indiv_vec[ii]
    if (!grepl("sex", sex_indiv)) {
      sex = rep(ctr, nrow(phen))
      ctr = ctr - 1
    } else {
      sex = true[[sex_indiv]]
    }
    res[, ii + 2] = (true[[lia_indiv]] >= thr[sex + 1]) + 0L
  }#first calculated status of the siblings, then replacing the values with the indicator of sibling status and number of siblings 
  res[, ncol(res) - 1] = ifelse(res[,ncol(res) - 2] + res[,ncol(res) - 1] > 0, 1, 0 ) 
  res[, ncol(res) - 2] = 2
  res[["GWAX"]] = as.numeric(rowSums(res[,c("CHILD_STATUS", "P1_STATUS", "P2_STATUS", "SIB_STATUS")]) > 0)
  newDist = gsub("\\.phen", "_prev020_sex\\.phen", phen.files[i])
  fwrite(res, newDist, quote = FALSE, sep = " ")
}


Klist = list(
  matrix(prev,ncol = 2),
  outer(c(5,1)/5, prev, "*"),
  outer(5:1/5, prev, "*")
  #  prev * 10:1/10
)



for (i in seq_along(true.files))  {
  print(i)
  (thr <- qnorm(1 - Klist[[2]]))
  true = fread(true.files[i])

  input <- tibble(
    FID           = 1:nrow(true),
    IID           = 1:nrow(true),
    child_gen     = true$offspringgeno_lia,
    child_full    = true$offspring_lia,
    child_sex     = true$sex + 1,
    father_full   = true$parents_lia.1,
    mother_full   = true$parents_lia.2
  )
  # skipping siblings for now

  input[["father_status"]] <- assign_status(input$father_full, thr[,1])
  input[["mother_status"]] <- assign_status(input$mother_full, thr[,2])
  input[["child_status"]]  <- assign_status_sex(input$child_full, thr, input_sex = input$child_sex)
  
  omni_dat = omniscience(
    input = input,
    h2 = .5,
    thr = thr,
    NA_frac = .2,
    init_samp = 2e6,
    account_for_sex = TRUE
  )
  omni_dat = omni_dat[,-grep("string", colnames(omni_dat))]
  
  newDist = gsub("sibs2", "omni_sex_5grps_prev010_sibs2", phen.files[i])
  fwrite(omni_dat, newDist, sep = " ", quote = FALSE, na = "NA")
}




#phen.files = str_subset(phen.files, pattern = "/sibs")
#phen.files = str_subset(phen.files, pattern = "\\\\sibs")



ltfh = fread("D:\\Work\\Project1\\simulatedData\\LTFHsibs2_100kx10k_NF_C1000_v5_maf005_prev020_sex.phen")
#omni = fread("D:\\Work\\Project1\\simulatedData\\omni_sex_5grps_sibs2_100kx10k_NF_C1000_v1_maf005.phen")
omni = fread("D:\\Work\\Project1\\simulatedData\\omni_sex_5grps_sibs2_100kx10k_NF_C1000_v5_maf005.phen")
true = fread("D:\\Work\\Project1\\simulatedData\\sibs2_100kx10k_NF_C1000_v5_maf005.true")
cor(ltfh$ltfh_nosib, true$offspringgeno_lia)
da = left_join(x = omni, y = true, by = c("FID", "IID"))
cor(da$post_mean_liab, da$offspringgeno_lia)


#ph = fread("D:\\Work\\Project1\\simulatedData\\sibs2_100kx10k_NF_C1000_v9_maf005_prev020_sex.phen")
ph = fread("D:\\Work\\Project1\\simulatedData\\LTFHsibs2_100kx10k_NF_C1000_v9_maf005_prev020_sex.phen")
table(ph$ltfh_cc, sex = true$sex)
