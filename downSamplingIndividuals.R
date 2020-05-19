library(data.table)
library(stringr)

#get list of all the summary stats to use, and their true values.
#fileRoot = "/home/emp/faststorage/project1/simulatedData"
fileRoot = "D:\\Work\\Project1\\SimulatedData"
base = "sibs2"
sword = "100kx10k"
cword = "C500"
bword = "_NF"
identifier = "maf05f"

true.files = list.files(path = fileRoot, pattern = "\\.true", full.names = T)
true.files = str_subset(true.files, pattern = base)
true.files = str_subset(true.files, pattern = sword)
true.files = str_subset(true.files, pattern = cword)
true.files = str_subset(true.files, pattern = bword)
true.files = str_subset(true.files, pattern = identifier)


keeplist.files = gsub("true", "keeplist", true.files)


phen.files = list.files(path = fileRoot, pattern = "\\.phen", full.names = T)
phen.files = str_subset(phen.files, pattern = base)
phen.files = str_subset(phen.files, pattern = sword)
phen.files = str_subset(phen.files, pattern = cword)
phen.files = str_subset(phen.files, pattern = bword)
phen.files = str_subset(phen.files, pattern = identifier)

phen.files = str_subset(phen.files, pattern = "/sibs")


ExcludeListGen = function(keeplist.files, ratio = .5, n_tot) {
  for (i in seq_along(keeplist.files)) {
    phen = fread(phen.files[i])
    FID = 1:nrow(phen)
    
    phen_col = "CHILD_STATUS"
    ctrls = sample(FID[phen[[phen_col]] == 0], size = n_tot*ratio)
    cases = sample(FID[phen[[phen_col]] == 1], size = n_tot*(1 - ratio))
    all_indivs = sort(c(cases, ctrls))
    removers = data.frame("FID" = all_indivs, "IID" = all_indivs)
    fwrite(removers, keeplist.files[i], sep = " ", quote = F)
  }
}

ExcludeListGen(keeplist.files, ratio = .5, n_tot = 10000)
