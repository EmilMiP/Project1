source("D:/Work/Project1/Code/Project1/omni.R")
setwd("D:/Work/Project1/")

library(data.table)
library(stringr)

identifier = "prev010_maf05f"
sword      = "100kx10k"


true.files = list.files("./simulatedData", full.names = TRUE, pattern = "true")
true.files = str_subset(true.files, pattern = identifier)
true.files = str_subset(true.files, pattern = sword)

phen.files = list.files("./simulatedData", full.names = TRUE, pattern = "phen")
phen.files = str_subset(phen.files, pattern = identifier)
phen.files = str_subset(phen.files, pattern = sword)


ltfh.files = str_subset(phen.files, pattern = "LTFH")
phen.files = phen.files[grep("/sibs", phen.files)]


nsib = 0

prev = c(.08, .02) * 2
#prev = .1

K = outer(5:1/5, prev, "*")

(thr = qnorm(1 - K))


for (k in seq_along(phen.files)) {
  print(k)
  phen = fread(phen.files[k])
  true = fread(true.files[k])
  
  input = tibble(
    FID = phen$FID,
    IID = phen$IID,
    child_status  = assign_status_sex(true$offspring_lia, thr, true$sex + 1),
    child_sex     = true$sex + 1,
    father_status = assign_status(true$parents_lia.1, thr[,1]),
    mother_status = assign_status(true$parents_lia.2, thr[,2])
  )
  if (nsib > 0) {
    sib_cols = str_subset(colnames(true), "sib")
    sib_sex_cols = str_subset(sib_cols, "sex")
    sib_lia_cols = str_subset(sib_cols, "g_lia")
    if (nsib > length(sib_lia_cols)) {
      cat("trying to use too many siblings! \n")
    }
    for (i in 1:nsib) {
      input[[paste("sib", i, "_status", sep = "")]] = assign_status_sex(true[[sib_lia_cols[i]]], thr, true[[sib_sex_cols[i]]] + 1)
      input[[paste("sib", i, "_sex", sep = "")]]    = true[[sib_sex_cols[i]]] + 1
    }
  }
  omni_dat = omniscience(input = input, 
                         h2 = 0.5, 
                         thr = thr,
                         nsib = nsib,
                         account_for_sex = TRUE,
                         init_samp = 1e6,
                         NA_frac = 1/5)

  omni_dat = omni_dat[,-grep("string", colnames(omni_dat))]
  
  out_dist = gsub("/sibs2", "/omni_sex_2grp_sibs2", phen.files[k])
  
  fwrite(omni_dat, out_dist, sep = " ", quote = FALSE, na = "NA")
}


#downsample files
ds.files = gsub("phen", "keeplist", phen.files)

res = matrix(NA, ncol = 3, nrow = length(phen.files))
colnames(res) <- c("omni", "ltfh", "gwas")
for (i in seq_along(phen.files)) {
  phen = fread(phen.files[i])
  true = fread(true.files[i])
  ltfh = fread(ltfh.files[i])
  print(round(100 * cov(true[,c("offspring_lia", "offspringgeno_lia", "parents_lia.1", "parents_lia.2")]), 2))
  omni = fread(gsub("/sibs2", "/omni_sex_2grp_sibs2", phen.files[i]))
  ds = fread(ds.files[i])

  ph   = ds %>% left_join(omni, by = c("FID", "IID")) %>% left_join(true, by = c("FID", "IID")) %>% left_join(ltfh, c("FID", "IID"))
  res[i, 1] = cor(ph$post_mean_liab, ph$offspringgeno_lia)
  res[i, 2] = cor(ph$ltfh_nosib, ph$offspringgeno_lia)
  res[i, 3] = cor(, ph$offspringgeno_lia)
}
res

res[,1]/res[,2]

par(mfrow = c(2,1))
plot(ph$post_mean_liab, ph$offspringgeno_lia, pch = 19)
abline(b = 1, a = 0)
plot(ph$ltfh_nosib, ph$offspringgeno_lia)
abline(b = 1, a = 0)


library(ggplot2)
ph   = omni_dat %>% left_join(true) %>% left_join(ltfh)
ggplot(ph) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(post_mean_liab, offspringgeno_lia, color = as.factor(child_sex)), alpha = .3) +
  geom_abline(col = "black") + 
  theme(legend.position = "top")

ggplot(ph) +
  bigstatsr::theme_bigstatsr() + 
  geom_point(aes(ltfh_nosib, offspringgeno_lia, color = as.factor(child_sex)), alpha = .3) +
  geom_abline(col = "black") + 
  theme(legend.position = "top")




round(100 * cov(true[,c("offspring_lia", "offspringgeno_lia", "parents_lia.1", "parents_lia.2", "parents_geno_lia.1", "parents_geno_lia.2")]), 2)
#, "sibling_lia.1", "sibling_geno_lia.1", "sibling_lia.2", "sibling_geno_lia.2"


par(mfrow = c(2,5))
lias = str_subset(colnames(true), "_lia")
for (i in seq_along(lias)) {
  cur_lia = lias[i]
  print(var(true[[cur_lia]]))
  qqnorm(true[[cur_lia]])
  abline(b = 1, a = 0)
}

